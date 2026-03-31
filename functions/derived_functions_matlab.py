import numpy as np
import gsw
import PyCO2SYS as pyco2
import xarray as xr
# import functions.O2sol as O2sol # this function has an error
# import functions.MLD_4grp as fl_mld # unused, used the (very similar) script mld_calcs_campbell now
import functions.mld_calcs_campbell as mld_cb
import os.path
import matplotlib.pyplot as plt
import pandas as pd
import matlab.engine


# ### LIAR/LIPHR wrapper
# ## Functions for carbonate system calculations (pH, TALK etc)
# MLD Calculations

def sigma0(salinity,temperature,lon,lat,pressure):
    SA = gsw.SA_from_SP(salinity,
                        pressure,
                        lon,
                        lat)

    CT = gsw.CT_from_t(SA,
                       temperature,
                       pressure)

    sigma = gsw.sigma0(SA,CT)
    
    return sigma, CT


def spiciness0(salinity,temperature,lon,lat,pressure):
    SA = gsw.SA_from_SP(salinity,
                        pressure,
                        lon,
                        lat)

    CT = gsw.CT_from_t(SA,
                       temperature,
                       pressure)

    spiciness = gsw.spiciness0(SA,CT)
    
    return spiciness

def calc_mld_wrapper(ds, suffix):
    mld_deboyer_interp=np.zeros(len(ds['N_PROF']))
    for k in range(len(ds.N_PROF.values)):
        par_subset = ['sigma0','PRES'+suffix,'depth']
        prof_data = ds[par_subset].isel(N_PROF=k).sortby('PRES'+suffix).dropna(subset=['sigma0'],dim='N_LEVELS')
        sigma_p = prof_data.sigma0.values
        dep_p = prof_data.depth.values

        if (len(dep_p)==0): # this is also the case when sigma0 is absent
            continue
        if (len(dep_p) > 0 and dep_p[0] > 30): # do not calculate MLD if shallowest depth is >30m (could be different value)
            continue 

        ref_depth = 10 if dep_p[0] <= 10 else 'shallowest'
        try:
            mld_deboyer_interp[k] = mld_cb.calc_mld(  
                                        dep_p,
                                        sigma_p,
                                        ref_depth=ref_depth,
                                        ref_reject=True,
                                        sigma_theta_crit=0.03,
                                        crit_method='interp',
                                        bottom_return='NaN',)
        except:
            mld_deboyer_interp[k] = np.nan
            pass
    ds['MLD'] = xr.DataArray(mld_deboyer_interp, coords={'N_PROF': ds['N_PROF']}, dims=['N_PROF'])
    return ds 


def LIPHR_matlab(LIPHR_path,Coordinates,Measurements,MeasIDVec,OAAdjustTF=False,VerboseTF=False):
#launch MATLAB engine API
    eng = matlab.engine.start_matlab()

    #convert inputs to MATLAB double
    Measurements = matlab.double([Measurements])
    Coordinates = matlab.double([Coordinates])
    MeasIDVec = matlab.double([MeasIDVec])
    
    #squeeze
    Measurements = eng.squeeze(Measurements)
    Coordinates = eng.squeeze(Coordinates)

    #need to make sure LIAR subfolders added to matlab path
    eng.addpath(eng.genpath(LIPHR_path))

    #call MATLAB function
    results = eng.LIPHR(Coordinates,Measurements,MeasIDVec,'OAAdjustTF', OAAdjustTF, 'VerboseTF', VerboseTF)
    eng.quit()

    results = np.asarray(results)   
    return results


def ESPER_matlab(LIPHR_path,DesiredVariables,Coordinates,Measurements,MeasIDVec_ESPER,Equations, Dates, ESPER_type, VerboseTF, pHCalcTF):
    # ESPER_type can be MX, NN, or LIR
    #launch MATLAB engine API
    eng = matlab.engine.start_matlab()

    #convert inputs to MATLAB double
    Measurements = matlab.double([Measurements])
    Coordinates = matlab.double([Coordinates])
    MeasIDVec_ESPER = matlab.double([MeasIDVec_ESPER])
    Equations = matlab.double([Equations])
    Dates = matlab.double([Dates])
    DesiredVariables = matlab.double([DesiredVariables])

    #squeeze
    Measurements = eng.squeeze(Measurements)
    Coordinates = eng.squeeze(Coordinates)

    #need to make sure LIAR subfolders added to matlab path
    eng.addpath(eng.genpath(LIPHR_path))

    #call MATLAB function
    results = eng.ESPER_wrapper_for_python(DesiredVariables,Coordinates,Measurements,MeasIDVec_ESPER,Equations,
                        Dates, ESPER_type, VerboseTF, pHCalcTF)
    eng.quit()

    results = np.asarray(results)   
    return results


def LIAR_matlab(LIAR_path,Coordinates,Measurements, MeasIDVec, pres_name, temp_name, sal_name,VerboseTF=False):
#launch MATLAB engine API
    eng = matlab.engine.start_matlab()

    #convert inputs to MATLAB double
    Measurements = matlab.double([Measurements])
    Coordinates = matlab.double([Coordinates])
    MeasIDVec = matlab.double([MeasIDVec])

    #squeeze
    Measurements = eng.squeeze(Measurements)
    Coordinates = eng.squeeze(Coordinates)
    
    #need to make sure LIAR subfolders added to matlab path
    eng.addpath(eng.genpath(LIAR_path))

    #call MATLAB function
    results = eng.LIAR(Coordinates,Measurements,MeasIDVec,'VerboseTF',VerboseTF)
    eng.quit()

    #convert matlab double output back to numpy array
    results = np.asarray(results)

    return results

def LINR_matlab(LIAR_path,Coordinates,Measurements, MeasIDVec,VerboseTF=False):
#launch MATLAB engine API
    eng = matlab.engine.start_matlab()

    #convert inputs to MATLAB double
    Measurements = matlab.double([Measurements])
    Coordinates = matlab.double([Coordinates])
    MeasIDVec = matlab.double([MeasIDVec])

    #squeeze
    Measurements = eng.squeeze(Measurements)
    Coordinates = eng.squeeze(Coordinates)
    
    #need to make sure LIR subfolders added to matlab path
    eng.addpath(eng.genpath(LIAR_path))

    #call MATLAB function
    results = eng.LINR(Coordinates,Measurements,MeasIDVec,'VerboseTF',VerboseTF)
    eng.quit()

    #convert matlab double output back to numpy array
    results = np.asarray(results)

    return results

def calculate_carbonate_parameters(argo_n, matlab_code_dir, pres_name, temp_name, sal_name, data_type_to_process, verbose):

    if verbose:
        print('using ' + data_type_to_process + ' variables')
    
    #initialise pH 25c and DIC variables
    if verbose:
        print('initializing empty variables')
    # argo_n['TALK_LIAR' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    # argo_n['TALK_LIAR' + data_type_to_process][:] = np.nan
    # calculating a number of different TA estimates for now to compare them all 
    argo_n['TALK_ESPER_NN' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['TALK_ESPER_NN' + data_type_to_process][:] = np.nan
    argo_n['TALK_ESPER_LIR' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['TALK_ESPER_LIR' + data_type_to_process][:] = np.nan
    argo_n['TALK_ESPER_MX' + data_type_to_process]= (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['TALK_ESPER_MX' + data_type_to_process][:] = np.nan

    argo_n['PH_25C_TOTAL' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['PH_25C_TOTAL' + data_type_to_process][:] = np.nan
    argo_n['PH_25C_0db_TOTAL' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['PH_25C_0db_TOTAL' + data_type_to_process][:] = np.nan

    # argo_n['DIC_LIAR' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    # argo_n['DIC_LIAR' + data_type_to_process][:] = np.nan
    # argo_n['pCO2_LIAR_W17' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    # argo_n['pCO2_LIAR_W17' + data_type_to_process][:] = np.nan
    
    argo_n['DIC_ESPER_MX' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['DIC_ESPER_MX' + data_type_to_process][:] = np.nan
    argo_n['PCO2_ESPER_MX'+ data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['PCO2_ESPER_MX' + data_type_to_process][:] = np.nan
 
    argo_n['DIC_ESPER_MX_THERMO' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['DIC_ESPER_MX_THERMO' + data_type_to_process][:] = np.nan
    argo_n['PCO2_ESPER_MX_THERMO'+ data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['PCO2_ESPER_MX_THERMO' + data_type_to_process][:] = np.nan
 
    argo_n['PCO2_ESPER_MX_W17'+ data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['PCO2_ESPER_MX_W17' + data_type_to_process][:] = np.nan
 
    argo_n['DIC_ESPER_NN' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['DIC_ESPER_NN' + data_type_to_process][:] = np.nan
    argo_n['PCO2_ESPER_NN_W17' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['PCO2_ESPER_NN_W17' + data_type_to_process][:] = np.nan

    argo_n['DIC_ESPER_LIR' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['DIC_ESPER_LIR' + data_type_to_process][:] = np.nan
    argo_n['PCO2_ESPER_LIR_W17' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['PCO2_ESPER_LIR_W17' + data_type_to_process][:] = np.nan

    argo_n['pH_insitu_corr' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    argo_n['pH_insitu_corr' + data_type_to_process][:] = np.nan
    # argo_n['pH_25C_corr'] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
    # argo_n.pH_25C_corr[:] = np.nan
    argo_n['bias_corr' + data_type_to_process] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof
    argo_n['bias_corr' + data_type_to_process][:] = np.nan

    nprof_n = argo_n.sizes['N_PROF']
    ##### Calc float TALK       
    #repeat lats, lons to match pressure shape
    lons_rep = np.tile(argo_n.LONGITUDE.values,(argo_n[pres_name].shape[1],1)).T
    lats_rep = np.tile(argo_n.LATITUDE.values,(argo_n[pres_name].shape[1],1)).T

    #set Si and PO4 inputs
    #if nitrate, then use redfield for Si and PO4?, otherwise set to 0 
    if verbose:
        print('estimating silicate and phosphate')
    if 'NITRATE' + data_type_to_process in argo_n.keys():
        SI = argo_n['NITRATE' + data_type_to_process]*2.5
        SI = SI.where(~np.isnan(SI), 0)
        PO4 = argo_n['NITRATE' + data_type_to_process]/16
        PO4 = PO4.where(~np.isnan(PO4),0)
        # Coordinates = np.stack((lons_rep.flatten(), 
        #                 lats_rep.flatten(), 
        #                 argo_n[pres_name].values.flatten()),
        #                 axis=1)
        # Measurements = np.stack((argo_n[sal_name].values.flatten(), 
        #                     argo_n[temp_name].values.flatten(), 
        #                     argo_n.NITRATE_ADJUSTED.values.flatten(), 
        #                     argo_n.DOXY_ADJUSTED.values.flatten()),
        #                     axis=1)
        # MeasIDVec = [1, 7, 3, 6]

    else:
        SI = np.zeros((argo_n['PH_IN_SITU_TOTAL' + data_type_to_process].shape))
        PO4 = np.zeros((argo_n['PH_IN_SITU_TOTAL' + data_type_to_process].shape))

    if verbose:
        print('constructing coordinates and measurements to pass to Alkalinity algorithms')
    Coordinates = np.stack((lons_rep.flatten(), 
                    lats_rep.flatten(), 
                    argo_n[pres_name].values.flatten()),
                    axis=1)
    Measurements = np.stack((argo_n[sal_name].values.flatten(), 
                        argo_n[temp_name].values.flatten(),
                        argo_n['DOXY' + data_type_to_process].values.flatten()),
                        axis=1)
    MeasIDVec = [1, 7, 6]                            


    # if verbose:
    #     print('LIAR TA')
    # results = LIAR_matlab(matlab_code_dir,
    #                                         Coordinates.tolist(),
    #                                         Measurements.tolist(),
    #                                         MeasIDVec,
    #                                         pres_name, temp_name, sal_name,
    #                                         VerboseTF=False)                                  

    # argo_n['TALK_LIAR' + data_type_to_process] = (['N_PROF','N_LEVELS'],
    #                         np.reshape(np.asarray(results),argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))

    # ESPER Calculations (not planning to keep all)
    MeasIDVec_ESPER = [1, 2, 6] # S, T, O2 - different numbering than v2 LIRs
    Equations_ESPER = 7 # for ESPER - asking to use equation w/ S, T, and O2 only 
    DesiredVariables_ESPER = [1]
    # calculate decimal_year for ESPER
    da = argo_n.JULD
    decimal_year = da.dt.year + (da.dt.dayofyear - 1 + (da.dt.hour * 3600 + da.dt.minute * 60 + da.dt.second) / 86400) / (365 + da.dt.is_leap_year)
    dates_rep = np.tile(decimal_year,(argo_n[pres_name].shape[1],1)).T

    # calculate ESPER mixed
    results_TALK_ESPER_MX = ESPER_matlab(matlab_code_dir,
                                                    DesiredVariables_ESPER,
                                                    Coordinates.tolist(),
                                                    Measurements.tolist(),
                                                    MeasIDVec_ESPER,
                                                    Equations_ESPER, 
                                                    dates_rep.flatten().tolist(), 
                                                    'MX',
                                                    0, 0)
    argo_n['TALK_ESPER_MX' + data_type_to_process] = (['N_PROF','N_LEVELS'],
                            np.reshape(np.asarray(results_TALK_ESPER_MX),argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))
    argo_n['TALK_BGCArgoPlus'] = argo_n['TALK_ESPER_MX' + data_type_to_process].copy()

    # calculate ESPER NN
    results_TALK_ESPER_NN = ESPER_matlab(matlab_code_dir,
                                                    DesiredVariables_ESPER,
                                                    Coordinates.tolist(),
                                                    Measurements.tolist(),
                                                    MeasIDVec_ESPER,
                                                    Equations_ESPER, 
                                                    dates_rep.flatten().tolist(), 
                                                    'NN',
                                                    0, 0)
    argo_n['TALK_ESPER_NN' + data_type_to_process] = (['N_PROF','N_LEVELS'],
                            np.reshape(np.asarray(results_TALK_ESPER_NN),argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))
    
     # calculate ESPER LIR
    results_TALK_ESPER_LIR = ESPER_matlab(matlab_code_dir,
                                                    DesiredVariables_ESPER,
                                                    Coordinates.tolist(),
                                                    Measurements.tolist(),
                                                    MeasIDVec_ESPER,
                                                    Equations_ESPER, 
                                                    dates_rep.flatten().tolist(), 
                                                    'LIR',
                                                    0, 0)
    argo_n['TALK_ESPER_LIR' + data_type_to_process] = (['N_PROF','N_LEVELS'],
                            np.reshape(np.asarray(results_TALK_ESPER_LIR),argo_n.PH_IN_SITU_TOTAL_ADJUSTED.shape))
    
    if verbose:
        print('calculating pH at a constant 25C')
    ##### Calculate float pH at 25C, DIC and apply bias corr
    results = pyco2.sys(
            par1=argo_n['TALK_ESPER_MX' + data_type_to_process], 
            par2=argo_n['PH_IN_SITU_TOTAL' + data_type_to_process],
            par1_type=1,
            par2_type=3,
            temperature=argo_n[temp_name], 
            pressure=argo_n[pres_name], 
            salinity=argo_n[sal_name], 
            temperature_out=25.,#*np.ones(argo_n[pres_name].shape), #fixed 25C temperature
            pressure_out= argo_n[pres_name], 
            total_silicate=SI,
            total_phosphate=PO4,
            opt_pH_scale = 1, #total
            opt_k_carbonic=10, #Lueker et al. 2000
            opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
            opt_total_borate=2, # Lee et al. 2010
            opt_k_fluoride=2, # Perez and Fraga 1987
            opt_buffers_mode=1,
    )
    argo_n['PH_25C_TOTAL' + data_type_to_process] = (['N_PROF','N_LEVELS'], results['pH_total_out']) # pH normalized to 25C. still w/ in situ pressure

    if verbose:
        print('calculating pH at 25C and 0db')
    results = pyco2.sys(
            par1=argo_n['TALK_ESPER_MX' + data_type_to_process], 
            par2=argo_n['PH_IN_SITU_TOTAL' + data_type_to_process],
            par1_type=1,
            par2_type=3,
            temperature=argo_n[temp_name], 
            pressure=argo_n[pres_name], 
            salinity=argo_n[sal_name], 
            temperature_out=25.,#*np.ones(argo_n[pres_name].shape), #fixed 25C temperature
            pressure_out= 0, 
            total_silicate=SI,
            total_phosphate=PO4,
            opt_pH_scale = 1, #total
            opt_k_carbonic=10, #Lueker et al. 2000
            opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
            opt_total_borate=2, # Lee et al. 2010
            opt_k_fluoride=2, # Perez and Fraga 1987
            opt_buffers_mode=1,
    )
    argo_n['PH_25C_0db_TOTAL' + data_type_to_process] = (['N_PROF','N_LEVELS'], results['pH_total_out']) # pH normalized to 25C and 0db 
    
    # if verbose:
    #     print('calculating DIC LIAR and pCO2 without Williams bias correction')
    # # calculate pCO2 and DIC, using in situ temperature and pressure 
    # results = pyco2.sys(
    #         par1=argo_n['TALK_LIAR' + data_type_to_process], 
    #         par2=argo_n['PH_IN_SITU_TOTAL' + data_type_to_process],
    #         par1_type=1,
    #         par2_type=3,
    #         temperature=argo_n[temp_name], 
    #         pressure=argo_n[pres_name], 
    #         salinity=argo_n[sal_name], 
    #         temperature_out=argo_n['cons_temp'],
    #         pressure_out=0, 
    #         total_silicate=SI,
    #         total_phosphate=PO4,
    #         opt_pH_scale = 1, #total
    #         opt_k_carbonic=10, #Lueker et al. 2000
    #         opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
    #         opt_total_borate=2, # Lee et al. 2010
    #         opt_k_fluoride=2, # Perez and Fraga 1987
    #         opt_buffers_mode=1,
    # )

    # argo_n['DIC_LIAR' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['dic']) 
    #  
    # if verbose:
    #     print('calculating DIC LIAR and pCO2 without Williams bias correction')
    # # calculate pCO2 and DIC, using in situ temperature and pressure 
    # results = pyco2.sys(
    #         par1=argo_n['TALK_ESPER_MX' + data_type_to_process], 
    #         par2=argo_n['PH_IN_SITU_TOTAL' + data_type_to_process],
    #         par1_type=1,
    #         par2_type=3,
    #         temperature=argo_n[temp_name], 
    #         pressure=argo_n[pres_name], 
    #         salinity=argo_n[sal_name], 
    #         temperature_out=argo_n['cons_temp'],
    #         pressure_out=0, 
    #         total_silicate=SI,
    #         total_phosphate=PO4,
    #         opt_pH_scale = 1, #total
    #         opt_k_carbonic=10, #Lueker et al. 2000
    #         opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
    #         opt_total_borate=2, # Lee et al. 2010
    #         opt_k_fluoride=2, # Perez and Fraga 1987
    #         opt_buffers_mode=1,
    # )

    # argo_n['DIC_ESPER_MX' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['dic']) 
    # argo_n['pCO2_ESPER_MX' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['pCO2_out']) 

    if verbose:
        print('calculating DIC and pCO2 from ESPER MX without the Williams 2017 bias correction and without the Johnson et al. thermodynamic correction')
    # calculate pCO2 and DIC, using in situ temperature and pressure 
    results = pyco2.sys(
            par1=argo_n['TALK_ESPER_MX' + data_type_to_process], 
            par2=argo_n['PH_IN_SITU_TOTAL' + data_type_to_process],
            par1_type=1,
            par2_type=3,
            temperature=argo_n[temp_name], 
            pressure=argo_n[pres_name], 
            salinity=argo_n[sal_name], 
            temperature_out=argo_n[temp_name],
            pressure_out=argo_n[pres_name], 
            total_silicate=SI,
            total_phosphate=PO4,
            opt_pH_scale = 1, #total
            opt_k_carbonic=10, #Lueker et al. 2000
            opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
            opt_total_borate=2, # Lee et al. 2010
            opt_k_fluoride=2, # Perez and Fraga 1987
            opt_buffers_mode=1,
    )

    argo_n['DIC_ESPER_MX' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['dic'])  
    argo_n['PCO2_ESPER_MX' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['pCO2_out']) 
    argo_n['DIC_BGCArgoPlus'] = argo_n['DIC_ESPER_MX' + data_type_to_process].copy()
    argo_n['PCO2_BGCArgoPlus'] = argo_n['PCO2_ESPER_MX' + data_type_to_process].copy() 

    if verbose:
        print('calculating DIC and pCO2 from ESPER MX without the Williams 2017 bias correction and WITH the Johnson et al. thermodynamic correction')
    # calculate pCO2 and DIC, using in situ temperature and pressure 
    pk1 = -np.log10(results['k_carbonic_1'])
    pk2 = -np.log10(results['k_carbonic_2'])
    pk1_adjusted = pk1 - 0.014
    pk2_adjusted = pk2 + 0.014
    k1_adjusted = 10**(-pk1_adjusted)
    k2_adjusted = 10**(-pk2_adjusted)
    results = pyco2.sys(
            par1=argo_n['TALK_ESPER_MX' + data_type_to_process], 
            par2=argo_n['PH_IN_SITU_TOTAL' + data_type_to_process],
            par1_type=1,
            par2_type=3,
            temperature=argo_n[temp_name], 
            pressure=argo_n[pres_name], 
            salinity=argo_n[sal_name], 
            temperature_out=argo_n[temp_name],
            pressure_out=argo_n[pres_name], 
            total_silicate=SI,
            total_phosphate=PO4,
            opt_pH_scale = 1, #total
            opt_k_carbonic=10, #Lueker et al. 2000
            opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
            opt_total_borate=2, # Lee et al. 2010
            opt_k_fluoride=2, # Perez and Fraga 1987
            opt_buffers_mode=1,
            k_carbonic_1=k1_adjusted,
            k_carbonic_2=k2_adjusted
    )

    argo_n['DIC_ESPER_MX_THERMO' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['dic'])  
    argo_n['PCO2_ESPER_MX_THERMO' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['pCO2_out']) 
    # argo_n['DIC_BGCArgoPlus'] = argo_n['DIC_ESPER_MX_THERMO' + data_type_to_process].copy()
    # argo_n['PCO2_BGCArgoPlus'] = argo_n['PCO2_ESPER_MX_THERMO' + data_type_to_process].copy()

    if verbose:
        print('Calculating DIC ESPER NN')
    # calculate pCO2 and DIC, using in situ temperature and pressure 
    results = pyco2.sys(
            par1=argo_n['TALK_ESPER_NN' + data_type_to_process], 
            par2=argo_n['PH_IN_SITU_TOTAL' + data_type_to_process],
            par1_type=1,
            par2_type=3,
            temperature=argo_n[temp_name], 
            pressure=argo_n[pres_name], 
            salinity=argo_n[sal_name], 
            temperature_out=argo_n[temp_name],
            pressure_out=argo_n[pres_name], 
            total_silicate=SI,
            total_phosphate=PO4,
            opt_pH_scale = 1, #total
            opt_k_carbonic=10, #Lueker et al. 2000
            opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
            opt_total_borate=2, # Lee et al. 2010
            opt_k_fluoride=2, # Perez and Fraga 1987
            opt_buffers_mode=1,
    )

    argo_n['DIC_ESPER_NN' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['dic'])  

    if verbose:
        print('calculating DIC ESPER LIR')
    # calculate pCO2 and DIC, using in situ temperature and pressure 
    results = pyco2.sys(
            par1=argo_n['TALK_ESPER_LIR' + data_type_to_process], 
            par2=argo_n['PH_IN_SITU_TOTAL' + data_type_to_process],
            par1_type=1,
            par2_type=3,
            temperature=argo_n[temp_name], 
            pressure=argo_n[pres_name], 
            salinity=argo_n[sal_name], 
            temperature_out=argo_n[temp_name],
            pressure_out=argo_n[pres_name], 
            total_silicate=SI,
            total_phosphate=PO4,
            opt_pH_scale = 1, #total
            opt_k_carbonic=10, #Lueker et al. 2000
            opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
            opt_total_borate=2, # Lee et al. 2010
            opt_k_fluoride=2, # Perez and Fraga 1987
            opt_buffers_mode=1,
    )

    argo_n['DIC_ESPER_LIR' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['dic'])  

    if verbose:
        print('calculating pH correction')
    for p in range(nprof_n):
    
        # skip a profile if pH is above 10.  There seem to be pH's above 10 that causing 
        if any(argo_n['PH_IN_SITU_TOTAL' + data_type_to_process][p,:]>10) or all(np.isnan(argo_n['PH_IN_SITU_TOTAL' + data_type_to_process][p,:])):
            'print pH out of range'
            continue

        #apply pH bias correction   
        #find if there are valid values between fixed p levels 
    
        #if there are valid pressure levels and valid pH between 1480-1520 db, 
        #calc bias correction only in this depth band, if not, calc correction between 970 and 1520
        # if any((argo_n[pres_name][p,:]>1480) & (argo_n[pres_name][p,:]<1520)):
        if (~np.isnan(argo_n['PH_25C_TOTAL' + data_type_to_process][p,(argo_n[pres_name][p,:]>1480) & (argo_n[pres_name][p,:]<1520)])).any():
            inds = (argo_n[pres_name][p,:]>1480) & (argo_n[pres_name][p,:]<1520)
            correction = -0.034529*argo_n['PH_25C_TOTAL' + data_type_to_process][p,inds]+0.26709
            # if verbose:
            #     print(pres_name + ' found between 1480 and 1520, correction = ' + str(correction.values))
                   
        else:
            inds = (argo_n[pres_name][p,:]>900) & (argo_n[pres_name][p,:]<1520)
            correction = -0.034529*argo_n['PH_25C_TOTAL' + data_type_to_process][p,inds]+0.26709
            # if verbose:
            #     print(pres_name + ' between 970 and 1520, correction = ' + str(correction.values))
           
        # print(correction)
        # print(np.isnan(correction))
        # print(np.sum(~np.isnan(correction)))
        if np.sum(~np.isnan(correction.values))>0: # only do the correction if a non-nan value is present in "correction" 
            argo_n['bias_corr' + data_type_to_process][p] = np.nanmean(correction)
            argo_n['pH_insitu_corr' + data_type_to_process][p,:] = argo_n['PH_IN_SITU_TOTAL' + data_type_to_process][p,:]+argo_n['bias_corr' + data_type_to_process][p]
            # argo_n.pH_25C_corr[p,:] = argo_n.PH_25C_TOTAL_ADJUSTED[p,:]+argo_n.bias_corr[p]
        else:
            if verbose:
                print('no non-nan values for correction')
    #call CO2sys again to get pCO2 with corrected PH- do we need to include this here?
    if verbose:
        print(' Calculating pCO2')
    # results = pyco2.sys(
    #         par1=argo_n['TALK_LIAR' + data_type_to_process], 
    #         par2=argo_n['pH_insitu_corr' + data_type_to_process],
    #         par1_type=1,
    #         par2_type=3,
    #         temperature=argo_n[temp_name], 
    #         pressure=argo_n[pres_name], 
    #         salinity=argo_n[sal_name], 
    #         temperature_out=argo_n['cons_temp'],
    #         pressure_out=0, #argo_n[pres_name],
    #         total_silicate=SI,
    #         total_phosphate=PO4,
    #         opt_pH_scale = 1, #total
    #         opt_k_carbonic=10, #Lueker et al. 2000
    #         opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
    #         opt_total_borate=2, # Lee et al. 2010
    #         opt_k_fluoride=2, # Perez and Fraga 1987
    #         opt_buffers_mode=1,
    #         )

    # argo_n['pCO2_LIAR_W17' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['pCO2_out'])  

    results = pyco2.sys(
            par1=argo_n['TALK_ESPER_MX' + data_type_to_process], 
            par2=argo_n['pH_insitu_corr' + data_type_to_process],
            par1_type=1,
            par2_type=3,
            temperature=argo_n[temp_name], 
            pressure=argo_n[pres_name], 
            salinity=argo_n[sal_name], 
            temperature_out=argo_n['cons_temp'],
            pressure_out=0, #argo_n[pres_name],
            total_silicate=SI,
            total_phosphate=PO4,
            opt_pH_scale = 1, #total
            opt_k_carbonic=10, #Lueker et al. 2000
            opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
            opt_total_borate=2, # Lee et al. 2010
            opt_k_fluoride=2, # Perez and Fraga 1987
            opt_buffers_mode=1,
            )

    argo_n['PCO2_ESPER_MX_W17' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['pCO2_out'])  

    results = pyco2.sys(
            par1=argo_n['TALK_ESPER_NN' + data_type_to_process], 
            par2=argo_n['pH_insitu_corr' + data_type_to_process],
            par1_type=1,
            par2_type=3,
            temperature=argo_n[temp_name], 
            pressure=argo_n[pres_name], 
            salinity=argo_n[sal_name], 
            temperature_out=argo_n['cons_temp'],
            pressure_out=0, #argo_n[pres_name],
            total_silicate=SI,
            total_phosphate=PO4,
            opt_pH_scale = 1, #total
            opt_k_carbonic=10, #Lueker et al. 2000
            opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
            opt_total_borate=2, # Lee et al. 2010
            opt_k_fluoride=2, # Perez and Fraga 1987
            opt_buffers_mode=1,
            )

    argo_n['PCO2_ESPER_NN_W17' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['pCO2_out'])  

    results = pyco2.sys(
            par1=argo_n['TALK_ESPER_LIR' + data_type_to_process], 
            par2=argo_n['pH_insitu_corr' + data_type_to_process],
            par1_type=1,
            par2_type=3,
            temperature=argo_n[temp_name], 
            pressure=argo_n[pres_name], 
            salinity=argo_n[sal_name], 
            temperature_out=argo_n['cons_temp'],
            pressure_out=0, #argo_n[pres_name],
            total_silicate=SI,
            total_phosphate=PO4,
            opt_pH_scale = 1, #total
            opt_k_carbonic=10, #Lueker et al. 2000
            opt_k_bisulfate=1, # Dickson 1990 (Note, matlab co2sys combines KSO4 with TB. option 3 = KSO4 of Dickson & TB of Lee 2010)
            opt_total_borate=2, # Lee et al. 2010
            opt_k_fluoride=2, # Perez and Fraga 1987
            opt_buffers_mode=1,
            )

    argo_n['PCO2_ESPER_LIR_W17' + data_type_to_process] = (['N_PROF','N_LEVELS'],results['pCO2_out'])  
 
    return argo_n

# ### Calculate neutral density
def neutral_density(sal, temp, press, lon, lat):

    eng = matlab.engine.start_matlab()
    eng.addpath(eng.genpath('~/Documents/MATLAB/eos_legacy_gamma_n'))

    s = sal.shape
    if len(s) > 1: # for 2d data like from Argo
        nden = np.ones(s)*np.nan
        for p in range(s[0]):
            press[p,:][press[p,:]<0] = np.nan # remove negative pressures if they exist

            
            psal_matlab = matlab.double(sal[p,:].astype('float'))
            temp_matlab = matlab.double(temp[p,:].astype('float'))
            pres_matlab = matlab.double(press[p,:].astype('float'))
            lon_matlab = matlab.double(lon[p,:].astype('float'))
            lat_matlab = matlab.double(lat[p,:].astype('float'))

            # try:
            ndenl = eng.eos80_legacy_gamma_n(psal_matlab[0],temp_matlab[0],pres_matlab[0],lon_matlab[0],lat_matlab[0])
            # except:
                # print('Removing negative pressures')
                # press[p,:][press[p,:]<0] = np.nan # remove negative pressures
                # pres_matlab = matlab.double(press[p,:].astype('float'))
                # ndenl = eng.eos80_legacy_gamma_n(psal_matlab[0],temp_matlab[0],pres_matlab[0],lon_matlab[0],lat_matlab[0])

            # print(np.squeeze(ndenl))
            # ndenl = np.asarray(ndenl[0])
            ndenl = np.asarray(np.squeeze(ndenl))
            nden[p,:] = ndenl

    else: # from 1d data like from glodap
        psal_matlab = matlab.double(sal.astype('float'))
        temp_matlab = matlab.double(temp.astype('float'))
        pres_matlab = matlab.double(press.astype('float'))
        lon_matlab = matlab.double(lon.astype('float'))
        lat_matlab = matlab.double(lat.astype('float'))

        nden = eng.eos80_legacy_gamma_n(psal_matlab[0],temp_matlab[0],pres_matlab[0],lon_matlab[0],lat_matlab[0])
        nden = np.asarray(nden[0])[0]

    eng.quit()

    return nden

def apply_outlier_detection(argo_n, output_dir, matched_file_list, verbose):
    # reviewer_symbols = ['xr', 'ob', '+k', '^m', 'sg']
    group_outlier_dir= output_dir+ '../outlier_file_collection/'
    # print(argo_file[:7])

    wmo_n = str(argo_n['WMO_ID'].values)

    # wmo = argo_file[:7]
    # find all outlier files that match the WMO, if any
    # matched_file_list = [outlier_file for outlier_file in outlier_files if np.logical_and(wmo_n in outlier_file, not os.path.isdir(group_outlier_output_dir+ outlier_file))]
    # print(matched_file_list)
    # plot_changes = True

    # if len(matched_file_list)> 0: # only apply outliers if any are found 
    # argo_n = xr.open_dataset(output_output_dir+ argo_file)
    # create a copy to use for plotting removed points
    argo_n_orig = argo_n.copy(deep=True)

    outliers_removed = False
    # open all outlier files that match the WMO
    for file_n in matched_file_list: 
        if verbose:
            print('Looking at contents of outlier file ' + file_n)
        # file_n = matched_file_list[o]
        # print(file_n)
        with open(group_outlier_dir+file_n) as csvfile:
            df_out = pd.read_csv(csvfile)
            if verbose: print('File read in')
            var = df_out['Variable'].values
            nprof = df_out['N_PROF'].values
            nlevel = df_out['N_LEVELS'].values
            if verbose: print(f'{len(nprof)} outliers found in file ' + file_n)
            if len(nprof) > 0: 
                df_out['reviewer_initials'] = file_n.split('_')[4]
                if 'df_all' in locals():
                    df_all = pd.concat([df_all, df_out])
                else:
                    df_all = df_out
                if verbose: print('Applying outlier removal from file: ' + file_n)
                outliers_removed=True # sets to True if any values are being changed
                for i in range(0, len(nprof)): # go through each row, setting all values to nans
                    # Replace the data with nans
                    if type(var[i])==str:
                        # if var[i] ends in '_RO', change 'RO' to 'BGCArgoPlus' - deals with naming convention change
                        if var[i].endswith('_RO'):
                            var_name_for_outlier = var[i][:-2] + 'BGCArgoPlus'
                        elif var[i].endswith('_BGCArgoPlus'):
                            var_name_for_outlier = var[i]
                        else: 
                            print('Skipping variable ' + var[i] + ' in file:  ' + file_n + ' - does not end with _RO or _BGCArgoPlus')
                            continue # skip this variable if it doesn't end with _RO or _BGCArgoPlus - don't want to change original adjusted data
                        if verbose: print(var_name_for_outlier)
                        
                        argo_n[var_name_for_outlier].loc[{'N_PROF':int(nprof[i]), 'N_LEVELS':int(nlevel[i])}] = np.nan
                        
                        # if verbose: 
                        #     print('Adding OR to outlier flag for ' + var_name_for_outlier +'_flag: ')
                        #     print(argo_n[var_name_for_outlier +'_flag'].values.item())
        # if an outlier file exists for a float, that indicates that outlier detection was performed even if no outliers were removed
        # Add 'OR' to BGC variables even if no outliers were removed
        for var in ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOXY_ADJUSTED', 'NITRATE_ADJUSTED', 'PH_IN_SITU_TOTAL_ADJUSTED']:
            if var in argo_n.variables:
                var_name_for_outlier = var + '_BGCArgoPlus'
                if 'OR' not in argo_n[var_name_for_outlier +'_flag'].values.item():
                    if verbose: 
                        print('Adding OR to outlier flag for ' + var_name_for_outlier +'_flag: ')
                    argo_n[var_name_for_outlier +'_flag'] = argo_n[var_name_for_outlier +'_flag'].values.item() + 'OR_'
                    if verbose:
                        print(argo_n[var_name_for_outlier +'_flag'].values.item())

    # # if there were outliers removed and plot_changes is set to True, then plot profiles that were removed     
    # if np.logical_and(outliers_removed, plot_changes):

    #     # make a removed profile directory for this float if one doesn't exist
    #     float_RO_profile_dir= group_outlier_dir+ '/' + wmo_n + '/'
    #     if not os.path.isdir(float_RO_profile_dir):
    #         os.mkdir(float_RO_profile_dir)

    #     outlier_profiles = np.sort(df_all['N_PROF'].unique())
    #     for idx_p in range(0, len(outlier_profiles)):
    #         if np.isnan(outlier_profiles[idx_p]):
    #             continue 
    #         plot_filename = wmo_n + '_profile_' +  str(outlier_profiles[idx_p])

    #         prof_index = df_all['N_PROF']==outlier_profiles[idx_p]
    #         var_all = df_all['Variable'][prof_index]
    #         unique_variables = np.unique(var_all).tolist()
    #         if 'TEMP_ADJUSTED_RO' not in unique_variables:
    #             unique_variables.append('TEMP_ADJUSTED_RO')
    #         if 'PSAL_ADJUSTED_RO' not in unique_variables:
    #             unique_variables.append('PSAL_ADJUSTED_RO')                
    #         unique_reviewers = np.unique(df_all['reviewer_initials'])

    #         fig = plt.figure(figsize=(len(unique_variables)*10,10),)

    #         for idx_v in range(0, len(unique_variables)):


    #             # loop through all removed variables types for that profile
    #             # first plot that profile, plus ones before and after:
    #             ax = fig.add_subplot(1, len(unique_variables), idx_v+1)

    #             # make sure we don't exceed the limits of available profiles
    #             first_profile = int(outlier_profiles[idx_p]-2)
    #             last_profile = int(outlier_profiles[idx_p]+3)
    #             if first_profile<0:
    #                 first_profile=0
    #             if last_profile>len(argo_n_orig['N_PROF'].values)-1:
    #                 last_profile = argo_n_orig['N_PROF'][-1].values

    #             for prof_plot in range(first_profile,last_profile):
    #                 # print(prof_plot)
    #                 ax.plot(argo_n_orig[unique_variables[idx_v]].loc[{'N_PROF':prof_plot}].values, \
    #                     argo_n_orig['PRES'].loc[{'N_PROF':prof_plot}].values, label = prof_plot)
                    
    #             # plt.plot(argo_n_orig[unique_variables[idx_v]].loc[{'N_PROF':outlier_profiles[idx_p]}].values, \
    #             #          argo_n_orig['PRES'].loc[{'N_PROF':outlier_profiles[idx_p]}].values, label =)
    #             ax.set_title(unique_variables[idx_v] + ' Profile: ' + str(outlier_profiles[idx_p]))

    #             # plot removed variables:

    #             df_prof_var = df_all[prof_index].where(df_all['Variable'][prof_index]==unique_variables[idx_v]) 
    #             df_prof_var = df_prof_var.drop_duplicates().dropna(how='all').reset_index()

    #             # add an index item for the reviewers present in the profile
    #             prof_reviewers = np.unique(df_prof_var['reviewer_initials'])
    #             for pr in range(0, len(prof_reviewers)):
    #                 reviewer_match = unique_reviewers==prof_reviewers[pr]
    #                 rev_index = [i for i, x in enumerate(reviewer_match) if x] # get an index for the reviewer so that you have different legend values 

    #                 ax.plot(np.nan, np.nan, reviewer_symbols[rev_index[0]], label=unique_reviewers[rev_index])
    #             label=df_prof_var['reviewer_initials']
    #             # loop through rows and plot outlier points
    #             for RO_n in range(0, len(df_prof_var)):
    #                 reviewer_match = unique_reviewers==df_prof_var['reviewer_initials'][RO_n] 
    #                 rev_index = [i for i, x in enumerate(reviewer_match) if x] # get an index for the reviewer so that you have different legend values 
                    
    #                 # print(argo_n_orig[unique_variables[idx_v]].loc[{'N_PROF':outlier_profiles[idx_p], 'N_LEVELS':int(df_prof_var['N_LEVELS'][RO_n])}].values)
    #                 ax.plot(argo_n_orig[unique_variables[idx_v]].loc[{'N_PROF':int(outlier_profiles[idx_p]), 'N_LEVELS':int(df_prof_var['N_LEVELS'][RO_n])}].values, \
    #                     argo_n_orig['PRES'].loc[{'N_PROF':int(outlier_profiles[idx_p]), 'N_LEVELS':int(df_prof_var['N_LEVELS'][RO_n])}].values, reviewer_symbols[rev_index[0]] )
                    
    #             ax.invert_yaxis()
    #             ax.legend()
    #         plt.tight_layout()
    #         plt.savefig(f'{float_RO_profile_dir}{plot_filename}.png')
    #         plt.close(fig)
    argo_n_orig.close()
    # # save intermediate file
    # print('here')
    # argo_n.to_netcdf(output_output_dir+ argo_file[:-3] + '_RO.nc')
    # argo_n.close()        
    
    return argo_n

def calculate_derived_parameters(output_dir, file, matlab_code_dir, data_type_to_process, outlier_list, outlier_df, verbose=False):
# pass interpolated dataset and non-interpolated dataset through calculations of: 
#sigma0, conservative temperature (potential temperature) spiciness0, 
# dir_types = [output_dir, interpolated_directory]
# dir_types = [output_dir]
    try:
        pres_name = 'PRES' + data_type_to_process
        temp_name = 'TEMP' + data_type_to_process
        sal_name = 'PSAL' + data_type_to_process
        print(file)

        if verbose:
            print('opening "filtered" file')
        argo_n = xr.open_dataset(output_dir+ file)
        if verbose:
            print(file + ' opened')
        if 'WMO_ID' not in argo_n.keys():
            return
        wmo_n = argo_n['WMO_ID'].values
        if np.sum(~np.isnan(argo_n[temp_name]))==0:
            temp_name = 'TEMP_ADJUSTED'
            if np.sum(~np.isnan(argo_n[temp_name]))==0:
                temp_name = 'TEMP'
        if np.sum(~np.isnan(argo_n[sal_name]))==0:
            sal_name = 'PSAL_ADJUSTED'
            if np.sum(~np.isnan(argo_n[sal_name]))==0:
                sal_name = 'PSAL'
        if np.sum(~np.isnan(argo_n[pres_name]))==0:
            pres_name = 'PRES_ADJUSTED'
            if np.sum(~np.isnan(argo_n[pres_name]))==0:
                pres_name = 'PRES'
                if np.sum(~np.isnan(argo_n[pres_name]))==0:
                    print('No valid "' + pres_name + '" for ' + str(wmo_n) + ', skipping creation of a processed file.')
                    if os.path.isfile(output_dir+str(wmo_n)+'_Sprof_BGCArgoPlus_full.nc'): # delete file if it already exists
                        print('Existing processed file found, deleting')
                        os.remove(output_dir+str(wmo_n)+'_Sprof_BGCArgoPlus_full.nc')
                    return
        

        # get list of outliers and remove
        if verbose:
            print('getting list of outlier files')
        matched_df = outlier_df.where(outlier_df['wmo']==wmo_n).dropna(how='all')
        if verbose:
            print('Matched df:')
            print(matched_df)
        if len(matched_df)>0:
            matched_file_list = outlier_list[matched_df.index]
            if verbose:
                print('loading and removing outliers from csv files')
            argo_n = apply_outlier_detection(argo_n, output_dir, matched_file_list, verbose)
        else:
            if verbose:
                print('No outlier files found, skipping removal')
            

        if verbose:
            print('Calculating potential density, spiciness, conservative temperature')
        argo_n['sigma0'] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
        argo_n.sigma0[:] = np.nan

        argo_n['spiciness0'] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
        argo_n.spiciness0[:] = np.nan
        # start_time = time.perf_counter()

        lat_grid = np.repeat(argo_n.LATITUDE.values[:, np.newaxis], len(argo_n.N_LEVELS), axis=1)
        lon_grid = np.repeat(argo_n.LONGITUDE.values[:, np.newaxis], len(argo_n.N_LEVELS), axis=1)

        sigma0_temp, cons_Temp_temp = sigma0(argo_n[sal_name].values, 
                        argo_n[temp_name].values,
                        lon_grid,
                        lat_grid,
                        argo_n[pres_name].values)

        sigma0_da = xr.DataArray(sigma0_temp, dims=('N_PROF', 'N_LEVELS'), coords={'N_LEVELS': argo_n.N_LEVELS, 'N_PROF': argo_n.N_PROF})
        cons_temp_da = xr.DataArray(cons_Temp_temp, dims=('N_PROF', 'N_LEVELS'), coords={'N_LEVELS': argo_n.N_LEVELS, 'N_PROF': argo_n.N_PROF})

        # Add the potential density and conservative temperature to the Dataset
        argo_n['sigma0'] = sigma0_da
        argo_n['cons_temp'] = cons_temp_da

        # calculate spiciness
        spiciness0_temp = spiciness0(argo_n[sal_name].values, 
                        argo_n[temp_name].values,
                        lon_grid,
                        lat_grid,
                        argo_n[pres_name].values)

        spiciness_da = xr.DataArray(spiciness0_temp, dims=('N_PROF', 'N_LEVELS'), coords={'N_LEVELS': argo_n.N_LEVELS, 'N_PROF': argo_n.N_PROF})

        # Add spiciness to the Dataset
        argo_n['spiciness0'] = spiciness_da

        if verbose:
            print('Calculating neutral density')
        # crate arrays for lon and lat to pass into neutral density function 
        s = argo_n[sal_name].shape[1]
        lons_array = np.array([[x] * s for x in argo_n.LONGITUDE.values])
        lats_array = np.array([[x] * s for x in argo_n.LATITUDE.values])

        # neutral density
        nden_temp = neutral_density(argo_n[sal_name].values, argo_n[temp_name].values, argo_n[pres_name].values, lons_array, lats_array)
        nden_temp_da = xr.DataArray(nden_temp, dims=('N_PROF', 'N_LEVELS'), coords={'N_LEVELS': argo_n.N_LEVELS, 'N_PROF': argo_n.N_PROF})

        # # Add neutral density to the Dataset
        argo_n['gamma'] = nden_temp_da

        if verbose:
            print('Calculating depth')
        # Calculate depth from pressure (convert negative depth to positive)
        depth_temp = gsw.conversions.z_from_p(argo_n[pres_name].values,lat_grid)*(-1)
        depth_da = xr.DataArray(depth_temp, dims=('N_PROF', 'N_LEVELS'), coords={'N_LEVELS': argo_n.N_LEVELS, 'N_PROF': argo_n.N_PROF})
        argo_n['depth'] = depth_da

        if verbose:
            print('Calculating MLD')
        # calculate mld
        argo_n = calc_mld_wrapper(argo_n, suffix=data_type_to_process)
        
        # calculate carbonate system parameters
        if 'PH_IN_SITU_TOTAL_ADJUSTED' in argo_n.keys() and np.any(~np.isnan(argo_n.PH_IN_SITU_TOTAL_ADJUSTED)):
            if verbose:     
                print('Calculating carbonate system parameters')
            argo_n = calculate_carbonate_parameters(argo_n, matlab_code_dir, pres_name, temp_name, sal_name, data_type_to_process, verbose)
            # # temporarily running a second time, using the "Adjusted" data instead of the "ADJUSTED_RO" data
            # # commenting out for now - would have to rerun flags without the other automated outlier detection things for actual comparison. Can implement once I've tested the "_RO" run
            # argo_n = calculate_carbonate_parameters(argo_n, matlab_code_dir, pres_name, temp_name, sal_name, '_ADJUSTED', verbose)

        if verbose:
            print('Calculating O2 Sat. Conc.')
        argo_n['DOXY_SAT'] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
        argo_n.DOXY_SAT[:] = np.nan
        # argo_n['DOXY_SAT'] = O2sol.O2sol(argo_n[sal_name], argo_n[temp_name])
        argo_n['DOXY_SAT'] = gsw.O2sol_SP_pt(argo_n[sal_name], argo_n[temp_name])

        if verbose:
            print('Saving processed file')
        if os.path.isfile(output_dir+str(wmo_n)+'_Sprof_BGCArgoPlus_full.nc'): # delete file if it already exists
                os.remove(output_dir+str(wmo_n)+'_Sprof_BGCArgoPlus_full.nc')
        argo_n.to_netcdf(output_dir+str(wmo_n)+'_Sprof_BGCArgoPlus_full.nc')

        # save a trimmed version with many non-essential variables removed
        if verbose:
            print('Saving pared down processed file')
        strings_in_vars_to_drop = ['_QC_n_3_removed', '_QC_n_4_removed', '_profile_mode', 'ESPER', 'pH_insitu_corr', 'bias_corr']
        for string in strings_in_vars_to_drop:
            if verbose:
                print('Dropping variables with "' + string + '" in the name')
            vars_to_drop = [var for var in argo_n.variables if string in var]
            if len(vars_to_drop)>0:
                argo_n = argo_n.drop_vars(vars_to_drop)
                if verbose: print('Dropped variables with : ' + string) 
        if verbose:
            print('attempting to save trimmed file')
        if os.path.isfile(output_dir+ '../for_external_sharing/' + str(wmo_n)+'_Sprof_BGCArgoPlus.nc'): # delete file if it already exists
                if verbose: print('Old trimmed file exists, deleting')
                os.remove(output_dir+ '../for_external_sharing/' + str(wmo_n)+'_Sprof_BGCArgoPlus.nc')
                if verbose: print('file deleted')
        argo_n.to_netcdf(output_dir+ '../for_external_sharing/' + str(wmo_n)+'_Sprof_BGCArgoPlus.nc')


        argo_n.close()        
        
    except:
        print( 'Error processing: ' + file)
    return


def calculate_carbonate_parameters_only(output_dir, file, matlab_code_dir, data_type_to_process, verbose=False):
    # pass interpolated dataset and non-interpolated dataset through calculations of: 
    #sigma0, conservative temperature (potential temperature) spiciness0, 

    try:
        pres_name = 'PRES' + data_type_to_process
        temp_name = 'TEMP' + data_type_to_process
        sal_name = 'PSAL' + data_type_to_process
        print(file)

        if verbose:
            print('opening "filtered" file')
            try:
                argo_n = xr.open_dataset(output_dir+ file)
            except Exception as e:
                print(f"Error opening file {file}: {e}")
                return
        if verbose:
            print(file + ' opened')
        if 'WMO_ID' not in argo_n.keys():
            if verbose:
                print('WMO_ID not found in ' + file + ', skipping creation of a processed file.')
            return
        wmo_n = argo_n['WMO_ID'].values
        if np.sum(~np.isnan(argo_n[temp_name]))==0:
            temp_name = 'TEMP_ADJUSTED'
            if np.sum(~np.isnan(argo_n[temp_name]))==0:
                temp_name = 'TEMP'
        if np.sum(~np.isnan(argo_n[sal_name]))==0:
            sal_name = 'PSAL_ADJUSTED'
            if np.sum(~np.isnan(argo_n[sal_name]))==0:
                sal_name = 'PSAL'
        if np.sum(~np.isnan(argo_n[pres_name]))==0:
            pres_name = 'PRES_ADJUSTED'
            if np.sum(~np.isnan(argo_n[pres_name]))==0:
                pres_name = 'PRES'
                if np.sum(~np.isnan(argo_n[pres_name]))==0:
                    print('No valid "' + pres_name + '" for ' + str(wmo_n) + ', skipping creation of a processed file.')
                    if os.path.isfile(output_dir+str(wmo_n)+'_Sprof_BGCArgoPlus_flags_mode_only.nc'): # delete file if it already exists
                        print('Existing processed file found, deleting')
                        os.remove(output_dir+str(wmo_n)+'_Sprof_BGCArgoPlus_flags_mode_only.nc')
                    return
        
        if verbose:
            print('Calculating conservative temperature')
        # argo_n['sigma0'] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
        # argo_n.sigma0[:] = np.nan

        # argo_n['spiciness0'] = (['N_PROF','N_LEVELS'],np.empty(argo_n[pres_name].shape)) #nprof x nlevel
        # argo_n.spiciness0[:] = np.nan
        # start_time = time.perf_counter()

        lat_grid = np.repeat(argo_n.LATITUDE.values[:, np.newaxis], len(argo_n.N_LEVELS), axis=1)
        lon_grid = np.repeat(argo_n.LONGITUDE.values[:, np.newaxis], len(argo_n.N_LEVELS), axis=1)

        sigma0_temp, cons_Temp_temp = sigma0(argo_n[sal_name].values, 
                        argo_n[temp_name].values,
                        lon_grid,
                        lat_grid,
                        argo_n[pres_name].values)

        # sigma0_da = xr.DataArray(sigma0_temp, dims=('N_PROF', 'N_LEVELS'), coords={'N_LEVELS': argo_n.N_LEVELS, 'N_PROF': argo_n.N_PROF})
        cons_temp_da = xr.DataArray(cons_Temp_temp, dims=('N_PROF', 'N_LEVELS'), coords={'N_LEVELS': argo_n.N_LEVELS, 'N_PROF': argo_n.N_PROF})

        # Add the potential density and conservative temperature to the Dataset
        # argo_n['sigma0'] = sigma0_da
        argo_n['cons_temp'] = cons_temp_da

        # calculate carbonate system parameters
        if 'PH_IN_SITU_TOTAL_ADJUSTED' in argo_n.keys() and np.any(~np.isnan(argo_n.PH_IN_SITU_TOTAL_ADJUSTED)):
            if verbose:     
                print('Calculating carbonate system parameters')
            argo_n = calculate_carbonate_parameters(argo_n, matlab_code_dir, pres_name, temp_name, sal_name, data_type_to_process, verbose)

        if verbose:
            print('Saving processed file')
        if os.path.isfile(output_dir+str(wmo_n)+'_Sprof_BGCArgoPlus_flags_mode_only.nc'): # delete file if it already exists
                os.remove(output_dir+str(wmo_n)+'_Sprof_BGCArgoPlus_flags_mode_only.nc')
        argo_n.to_netcdf(output_dir+str(wmo_n)+'_Sprof_BGCArgoPlus_flags_mode_only.nc')

        argo_n.close()        
    except:
        print( 'Error processing: ' + file)
    return
# BGC Argo+ processing

## BGC Argo float processing outline (Bushinsky et al., dataset paper):

**"Float_Processing_v2026_04.ipynb":**

0. Download data if desired 
1. Read in Sprof files

**Call "flag_sensor_info.py"**

2. Apply appropriate flags and set all non-delayed mode data (A and R) to nan. Save in new "[var]_ADJUSTED_BGCArgoPlus" variables
	- BGCArgoPlus flag per variable 
		- “DMonly” - only delayed mode allowed
		- “F” - DAC QC flags applied (Bad data flagged as NaN)
		- "S" - Surface removal
		- "Inv" = Density Inversion 
		- “OR” - Outlier Removed
		<!-- - “O2Bias” - O2 Bias correction applied (Bushinsky et al., 2025, Nachod et al., in prep) -->
		<!-- - “Thermo” - Thermodynamic correction applied to pK1/pK2 (Johnson et al., in prep) -->
	- Go through automatic density inversion removal
 	- Automatically remove surface values that shouldn't be in the Sprof files
	- Detect and remove bottom oxygen hooks
4. Load meta data, note calibration type for oxygen sensor (air, not air), load sensor types
2. Option to look for outliers (**"outlier_removal_single_float.py"**)
 	- Go through profiles to identify outliers to be removed	
	- Save list of bad data as a csv file

**Call "derived_functions_matlab.py"**

3. Read in and apply text file of bad data. 
   	- Finds all outlier files and sets all outlier points in "[var]_ADJUSTED_BGCArgoPlus" to NaN.
	- Calculate gamma, potential density, potential temperature, spice, O2 sat
4. 	Calculate MLD
	- Save attribute info about what MLD was used
5. Calculate derived carbonate system parameters if pH is present

6. Save out new netcdf files for each float:
	- save out file with all intermediate steps in : ../processed/ "[WMO]_Sprof_BGCArgoPlus_full.nc"
	- remove extraneous variables that most people donʻt want
	- save out "[WMO]_Sprof_BGCArgoPlus.nc" in ../processed/for_external_sharing/ --> This is what will be shared for dataset paper


**"Float_gridding_v2026_04.ipynb"**

8. Gridding, **Call "float_data_gridding.py"**:
	- Read in ../processed/: "[WMO]_Sprof_BGCArgoPlus.nc"
	- save out monthly files ../processed/for_external_sharing/gridded/monthly/
	- concatenates monthly into merged files ../processed/for_external_sharing/gridded/

## Packages required to run our processing code
We originally made a .yml file available for creating an environment capable of running everything in this repository. However, I've found that it is much quicker and easier to install packages one by one or in groups for this repository, so instead of a .yml file here is a list of what you need to install. Let me know if I'm missing something.

Via conda-forge
 - numpy
 - pandas
 - xarray=2024.01.1
 - netCDF4
 - matplotlib
 - jupyter
 - python=3.12.*
 - scipy
 - cartopy
 - tqdm # for helpful progress bars
 - pip
 - dash
 - seaborn

 Via pip:
   - PyCO2SYS==1.8.3
   - gsw=3.4.*
   - dash-bootstrap-components
   - matlabengine==9.13.9 # this is the version for matlab R2022b. May be challenging to install - can also try going to "/[matlabroot]\extern\engines\python" and running "python setup.py install"


## Example code for analyzing float data
Examples of different types of scripts to use for reading in and working with these data can be found here: https://github.com/Hi-Cycles/BGC_Argo_Plus_Code_Repository 

Scripts currently available:
- Float_Glodap_Obs_Density.ipynb
- Float_file_exploration.ipynb
- NCP_nitrate_drawdown.ipynb
- Seasonal_Cycles.ipynb


# BGC Argo+ processing

## BGC Argo float processing outline (Bushinsky et al., dataset paper):


0. Download data if desired 
1. Read in Sprof files
2. Apply appropriate flags and set all non-delayed mode data (A and R) to nan. Save in new "[var]_ADJUSTED_BGCArgoPlus" variables
	- BGCArgoPlus flag per variable 
		- “DMonly” - only delayed mode allowed
		- “F” - DAC QC flags applied (Bad data flagged as NaN)
		- "S" - Surface removal
		- "Inv" = Density Inversion 
		- “OR” - Outlier Removed
		- “O2Bias” - O2 Bias correction applied (Bushinsky et al., 2025, Nachod et al., in prep)
		- “Thermo” - Thermodynamic correction applied to pK1/pK2 (Johnson et al., in prep)
	- Go through automatic density inversion removal
 	- Automatically remove surface values that shouldn't be in the Sprof files
	- Detect and remove bottom oxygen hooks
4. Load meta data, note calibration type for oxygen sensor (air, not air), load sensor types
2. Option to look for outliers
 	- Go through profiles to identify outliers to be removed	
	- Save list of bad data as a csv file
3. Read in and apply text file of bad data. 
   	1. Finds all outlier files and sets all datapoints in "[var]_ADJUSTED_BGCArgoPlus" to NaN.
	5. Calculate gamma, potential density, potential temperature, spice, O2 sat
4. 	Calculate MLD
	1. **Save attribute info about what MLD was used**
5. Calculate derived carbonate system parameters if pH is present

6. Save out new netcdf files for each float:
	- save out file with all intermediate steps in : ../processed/ "[WMO]_Sprof_BGCArgoPlus_full.nc"
	- remove extraneous variables that most people donʻt want
	- save out "[WMO]_Sprof_BGCArgoPlus.nc" in ../processed/for_external_sharing/ --> This is what will be shared for dataset paper
8. Gridding:
	- Read in ../processed/: "[WMO]_Sprof_BGCArgoPlus.nc"
	- save out monthly files ../processed/for_external_sharing/gridded/monthly/
	- concatenates monthly into merged files ../processed/for_external_sharing/gridded/

## Example code for analyzing float data
Examples of different types of scripts to use for reading in and working with these data can be found here: https://github.com/Hi-Cycles/BGC_Argo_Plus_Code_Repository 

Scripts currently available:
- Float_Glodap_Obs_Density.ipynb
- Float_file_exploration.ipynb
- NCP_nitrate_drawdown.ipynb
- Seasonal_Cycles.ipynb

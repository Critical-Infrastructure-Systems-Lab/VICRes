**vic_inp_prep.py** --- VIC Input Preparation Python Script\
writen by Vu Trung Dung (dtvu2205@gmail.com)
________________________________________________________________________________

# Overview

This Python script extracts climate data (e.g., precipitation, maximum and minimum temperature, and wind speed) from NetCDF files and formats them as input for the Variable Infiltration Capacity (VIC) model. The script supports selecting a spatial domain, a range of years, and specific locations for output.
________________________________________________________________________________

# Requirements

Python 3.8+\
Python packages:\
os – for file and directory operations\
pandas – for working with tabular data\
xarray – for handling NetCDF climate datasets\
________________________________________________________________________________

# Directory Structure

Your project directory should look like:\
project_dir/\
├─ prec/        # Precipitation NetCDF files (e.g., chirps-v2.0.2022.days_p05.nc)\
├─ tmax/        # Maximum temperature NetCDF files (e.g., era5.02m.tmax.1200.22-23.nc)\
├─ tmin/        # Minimum temperature NetCDF files (e.g., era5.02m.tmin.0000.22-23.nc)\
├─ wind/        # Wind speed NetCDF files (e.g., era5.10m.wind.0600.22-23.nc)\
├─ loc_sample.txt # List of locations for VIC input (tab-separated lon/lat)\
├─ vic_inp/     # Output folder for VIC input files (auto-created by script)\
├─ vic_input_prep.py # Python script
________________________________________________________________________________

# Input Files

Precipitation:\
File example: chirps-v2.0.2022.days_p05.nc\
Variable name: precip

Maximum Temperature:\
File example: era5.02m.tmax.1200.22-23.nc\
Variable name: mx2t

Minimum Temperature:\
File example: era5.02m.tmin.0000.22-23.nc\
Variable name: mn2t

Wind Speed:\
File example: era5.10m.wind.0600.22-23.nc\
Variable name: fg10

Locations file (loc_sample.txt):\
Tab-separated file with columns:\
lon    lat\
90.123 15.456\
91.234 16.567
________________________________________________________________________________

# Configuration

At the top of the script, adjust the following parameters:

Spatial domain (latitude and longitude)\
min_lat, max_lat = 15, 35\
min_lon, max_lon = 90, 105

Years to process\
sta_year = 2022\
end_year = 2023

Directories\
prec_dir = "prec"\
tmax_dir = "tmax"\
tmin_dir = "tmin"\
wind_dir = "wind"\
outp_dir = "vic_inp"
________________________________________________________________________________

# How It Works

Data Preprocessing:

- Defines the spatial domain and year range.

- Reads all NetCDF files for precipitation and slices them to the selected domain.

- Reads ERA5 data for temperature and wind, inverting latitude if necessary.

VIC Input Extraction:

- Loops through each location in loc_sample.txt.

- Extracts nearest grid values for precip, tmax, tmin, and wind.

- Converts temperature from Kelvin to Celsius.

- Writes tab-separated files for VIC input in vic_inp/.

Export Output: For each location, a file named gf_<lon>_<lat> is created. Files contain 4 columns: prec, tmax, tmin, wind.

________________________________________________________________________________

# Notes

Ensure that NetCDF variable names match your dataset.

The script automatically creates the output folder if it does not exist.

Locations that do not exactly match the grid will use the nearest-neighbor selection.
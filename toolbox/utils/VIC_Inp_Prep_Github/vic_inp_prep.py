# LIBRARIES ___________________________________________________________________
import os          # For handling file paths
import pandas as pd # For handling tabular data (VIC input/output)
import xarray as xr # For working with NetCDF climate data

# PRE-PROCESSING _____________________________________________________________

# Define the spatial domain you want to process (latitude and longitude range)
min_lat, max_lat = 15, 35
min_lon, max_lon = 90, 105

# Define the range of years to process
sta_year = 2022
end_year = 2023

# Set directories for input data and output VIC files
main_dir = os.getcwd()
prec_dir = os.path.join(main_dir, "prec")   # Precipitation data folder
tmax_dir = os.path.join(main_dir, "tmax")   # Maximum temperature folder
tmin_dir = os.path.join(main_dir, "tmin")   # Minimum temperature folder
wind_dir = os.path.join(main_dir, "wind")   # Wind speed folder
outp_dir = os.path.join(main_dir, "vic_inp") # Output folder for VIC files

# Ensure the output directory exists
os.makedirs(outp_dir, exist_ok=True)

# PROCESSING __________________________________________________________________

# -------------------- Precipitation --------------------
print("Processing precipitation data: START")
prec = []  # List to store each year's data
years = range(sta_year, end_year+1)
for year in years:
    print("year: "+str(year))
    file_name = f"chirps-v2.0.{year}.days_p05.nc"  # File naming convention
    file_path = os.path.join(prec_dir, file_name)
    
    # Open NetCDF file with xarray
    data = xr.open_dataset(file_path)
    
    # Slice the data to the defined spatial domain
    sliced_data = data.sel({"latitude": slice(min_lat, max_lat),
                            "longitude": slice(min_lon, max_lon)})
    prec.append(sliced_data)

# Concatenate all years along the time dimension
prec = xr.concat(prec, dim="time")
print("Processing precipitation data: DONE")

# -------------------- Maximum Temperature --------------------
print("Processing maximum temperature data: START")
file_name = "era5.02m.tmax.1200.22-23.nc"
file_path = os.path.join(tmax_dir, file_name)

# Open NetCDF file
data = xr.open_dataset(file_path)

# ERA5 latitude may be from north to south, invert if needed
inverted_data = data.isel({"latitude": slice(None, None, -1)})

# Slice to defined spatial domain
tmax = inverted_data.sel({"latitude": slice(min_lat, max_lat),
                          "longitude": slice(min_lon, max_lon)})
print("Processing maximum temperature data: DONE")

# -------------------- Minimum Temperature --------------------
print("Processing minimum temperature data: START")
file_name = "era5.02m.tmin.0000.22-23.nc"
file_path = os.path.join(tmin_dir, file_name)

data = xr.open_dataset(file_path)
inverted_data = data.isel({"latitude": slice(None, None, -1)})
tmin = inverted_data.sel({"latitude": slice(min_lat, max_lat),
                          "longitude": slice(min_lon, max_lon)})
print("Processing minimum temperature data: DONE")

# -------------------- Wind Speed --------------------
print("Processing windspeed data: START")
file_name = "era5.10m.wind.0600.22-23.nc"
file_path = os.path.join(wind_dir, file_name)

data = xr.open_dataset(file_path)
inverted_data = data.isel({"latitude": slice(None, None, -1)})
wind = inverted_data.sel({"latitude": slice(min_lat, max_lat),
                          "longitude": slice(min_lon, max_lon)})
print("Processing windspeed data: DONE")

# EXTRACTING AND WRITING OUT DATA FOR VIC _____________________________________
print("Extracting and writing out data for VIC: START")

# Load locations (longitude, latitude) where VIC input files are needed
loc = pd.read_csv("loc_sample.txt", sep="\t", header=None, names=["lon", "lat"])

# Loop over each location
for index, row in loc.iterrows():
    lon = row['lon']
    lat = row['lat']
    print("lat: "+str(lat)+", lon: "+str(lon))    
    
    # Extract the nearest grid point data for each variable
    sel_prec = prec["precip"].sel(latitude=lat, longitude=lon, method="nearest")
    sel_tmax = tmax["mx2t"].sel(latitude=lat, longitude=lon, method="nearest")
    sel_tmin = tmin["mn2t"].sel(latitude=lat, longitude=lon, method="nearest")
    sel_wind = wind["fg10"].sel(latitude=lat, longitude=lon, method="nearest")
    
    # Create a pandas DataFrame with VIC-required variables
    df = pd.DataFrame({
        "prec": sel_prec.values,               # Precipitation (mm/day)
        "tmax": sel_tmax.values - 273.15,      # Max temperature converted from K to °C
        "tmin": sel_tmin.values - 273.15,      # Min temperature converted from K to °C
        "wind": sel_wind.values                # Wind speed (m/s)
    })

    # Define output filename and write to file (tab-separated, no header/index)
    filename = f"gf_{lon:.4f}_{lat:.4f}"  # Original format, no .txt extension
    file_path = os.path.join(outp_dir, filename)
    df.to_csv(file_path, sep="\t", index=False, header=False, float_format="%.6f")

print("Extracting and writing out data for VIC: DONE")

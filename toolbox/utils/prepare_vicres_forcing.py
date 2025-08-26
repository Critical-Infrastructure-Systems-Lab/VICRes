#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 10:51:59 2025

@author: Hisham Eldardiry
"""

import ee
import pandas as pd
from datetime import datetime, timedelta

# Authenticate and initialize Earth Engine
ee.Initialize()


# --- CONFIGURATION ---
lat, lon = 11.9189, 107.9069
start_date = '2022-01-01'
end_date = '2022-12-31'

point = ee.Geometry.Point([lon, lat])

# --- Generate daily list of dates ---
start_dt = datetime.strptime(start_date, "%Y-%m-%d")
end_dt = datetime.strptime(end_date, "%Y-%m-%d")
date_list = [(start_dt + timedelta(days=i)).strftime("%Y-%m-%d") 
             for i in range((end_dt - start_dt).days + 1)]


# --- CHIRPS Daily Precipitation ---
chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY") \
    .filterDate(start_date, end_date) \
    .filterBounds(point)

# --- ERA5-Land Hourly ---
era5 = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY") \
    .filterDate(start_date, end_date) \
    .filterBounds(point) \
    .select(['temperature_2m', 'u_component_of_wind_10m', 'v_component_of_wind_10m'])

# %% Add derived bands: temperature in C, wind speed
def process_hourly(img):
    temp_c = img.select('temperature_2m').subtract(273.15).rename('temp_C')
    u = img.select('u_component_of_wind_10m')
    v = img.select('v_component_of_wind_10m')
    wind = u.pow(2).add(v.pow(2)).sqrt().rename('wind_speed')
    return img.addBands([temp_c, wind])

era5 = era5.map(process_hourly)




# %% Function to create daily combined image
def daily_reduce(date_str):
    date = ee.Date(date_str)
    next_date = date.advance(1, 'day')

    # Handle missing CHIRPS
    chirps_filtered = chirps.filterDate(date, next_date)
    chirps_img = ee.Image(0).rename('precipitation')  # default image
    chirps_img = ee.Image(ee.Algorithms.If(
        chirps_filtered.size().gt(0),
        chirps_filtered.mean().select('precipitation'),
        chirps_img
    ))

    # ERA5: tmin, tmax, wind
    daily_imgs = era5.filterDate(date, next_date)
    tmin = daily_imgs.select('temp_C').min()
    tmax = daily_imgs.select('temp_C').max()
    wind = daily_imgs.select('wind_speed').mean()

    combined = chirps_img.addBands([tmax.rename('tmax'), tmin.rename('tmin'), wind.rename('wind')])
    return combined.set('date', date.format('YYYY-MM-dd'))

# Create image collection for daily combined records
daily_images = ee.ImageCollection([daily_reduce(d) for d in date_list])

# %% Extract data as list of features
def extract_record(img):
    values = img.reduceRegion(
        reducer=ee.Reducer.first(),
        geometry=point,
        scale=5000
    )
    return ee.Feature(None, values).set('date', img.get('date'))

feature_collection = ee.FeatureCollection(daily_images.map(extract_record))
features = feature_collection.getInfo()['features']


# %% Convert to DataFrame
df = pd.DataFrame([
    {
        'date': f['properties']['date'],
        'precip': f['properties'].get('precipitation'),
        'tmax': f['properties'].get('tmax'),
        'tmin': f['properties'].get('tmin'),
        'wind': f['properties'].get('wind')
    }
    for f in features
])

# %% Save Output
df[['precip', 'tmax', 'tmin', 'wind']].to_csv(
    f"data_{lat:.4f}_{lon:.4f}", sep=' ', index=False, header=False
)

print(f"Saved forcing file for lat {lat}, lon {lon}")

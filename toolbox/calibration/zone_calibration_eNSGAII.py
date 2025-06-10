#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  4 11:27:30 2025

@author: Hisham Eldardiry
"""

"""
Hydrological Model Calibration Script - Refactored MOdular Version
"""

# Sequential Zone-Based Calibration (Sequential Calibration)
# ----------------------------------------------------------
# This calibration strategy performs optimization zone by zone in a sequential manner.
# Key features:
# - Runs one optimization problem per zone (e.g., subbasin or gauged reach).
# - Starts from the most upstream zone and proceeds downstream.
# - Each downstream zone uses optimized parameters from the upstream zone.
# - Allows for localized tuning and improved accuracy at individual sites.
# - Useful for large, heterogeneous basins with spatial variability.

# %% ====================== LIBRARIES ======================
import os
import time
import shutil
import datetime
from datetime import date, timedelta
import numpy as np
import pandas as pd
import multiprocessing
from platypus import Real, Problem, EpsNSGAII, ProcessPoolEvaluator, Hypervolume, nondominated

# %% ====================== DEFINE USER DIRECTORIES ======================
# Base directory structure (modify these paths as needed)
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

# Main directories
CALIB_DIR = os.path.join(BASE_DIR, 'toolbox','calibration')
RUNOFF_DIR = os.path.join(BASE_DIR, 'runoff')
ROUTING_DIR = os.path.join(BASE_DIR, 'routing')

# Subdirectories
TRACKING_DIRS = {
    'performance': os.path.join(CALIB_DIR, 'performance_tracking'),
    'station': os.path.join(CALIB_DIR, 'station_tracking'), 
    'soil': os.path.join(CALIB_DIR, 'soil_tracking')
}

RUNOFF_SUBDIRS = {
    'parameters': os.path.join(RUNOFF_DIR, 'parameters'),
    'output': os.path.join(RUNOFF_DIR, 'output'),
    'model': os.path.join(RUNOFF_DIR, 'model')
}

ROUTING_SUBDIRS = {
    'parameters': os.path.join(ROUTING_DIR, 'parameters'),
    'output': os.path.join(ROUTING_DIR, 'output'),
    'model': os.path.join(ROUTING_DIR, 'model'),
    'obs_discharge': os.path.join(ROUTING_DIR, 'parameters', 'obs_discharge_1996_2005')
}

# Important files
CALIB_SETTINGS = os.path.join(CALIB_DIR, 'calibration_setup.txt')
ZONE_FILE = os.path.join(CALIB_DIR, 'zone.txt')
SOIL_FILE = os.path.join(RUNOFF_SUBDIRS['parameters'], 'soil.txt')
ASCII_TEMP = os.path.join(ROUTING_SUBDIRS['parameters'], 'flowdirection.txt')
ROUTING_CONFIG = os.path.join(ROUTING_SUBDIRS['parameters'], 'configuration.txt')

# %% ====================== HELPER FUNCTIONS ======================
def create_directories():
    """Create all required directories if they don't exist"""
    for dir_path in TRACKING_DIRS.values():
        os.makedirs(dir_path, exist_ok=True)

def validate_environment():
    """Check all system dependencies exist"""
    required_dirs = [RUNOFF_SUBDIRS['model'], ROUTING_SUBDIRS['model']]
    for d in required_dirs:
        if not os.path.exists(d):
            raise FileNotFoundError(f"Critical directory missing: {d}")
            
def write_zone_param_ascii(zone_df, ascii_temp_path, output_path, value_column):
    """Write zone parameters to ASCII file format"""
    with open(ascii_temp_path, 'r') as f:
        header_lines = [next(f) for _ in range(6)]
        header = {line.split()[0]: float(line.split()[1]) for line in header_lines}
        ascii_data = np.loadtxt(f)

    nrows = int(header['nrows'])
    ncols = int(header['ncols'])
    xllcorner = header['xllcorner']
    yllcorner = header['yllcorner']
    cellsize = header['cellsize']
    nodata = int(header['NODATA_value'])

    param_grid = np.full((nrows, ncols), nodata, dtype=np.float32)

    for _, row in zone_df.iterrows():
        lat = row['Col3']
        lon = row['Col4']
        val = row[value_column]

        col = int((lon - xllcorner) / cellsize)
        row_ = int((yllcorner + nrows * cellsize - lat) / cellsize)

        if 0 <= row_ < nrows and 0 <= col < ncols and ascii_data[row_, col] != nodata:
            param_grid[row_, col] = val

    with open(output_path, 'w') as f:
        for line in header_lines:
            f.write(line)
        np.savetxt(f, param_grid, fmt="%.4f")


def get_pop_ind(vic_count, generation, pop_size):
    """Calculate population individual number"""
    return ((vic_count) - (generation-1)*pop_size)

# %% ====================== METRICS CALCULATION ======================
def calculate_nse(obs, sim, handle_negatives=True):
    """
    Calculate Nash-Sutcliffe Efficiency (NSE) with optional negative flow handling.
    
    Args:
        obs: Array of observed values
        sim: 2D array of simulated values [..., [year, month, day, flow]]
        handle_negatives: If True, uses abs() on simulated flows (original behavior)
                         If False, uses standard NSE calculation
    
    Returns:
        NSE value (float)
    """
    obs_mean = np.mean(obs)
    sim_flows = np.abs(sim) if handle_negatives else sim
    
    numerator = np.sum((sim_flows - obs) ** 2)
    denominator = np.sum((obs - obs_mean) ** 2)
    
    return 1 - numerator / denominator

def calculate_trmse(obs, sim, handle_negatives=True):
    """
    Calculate Transformed Root Mean Squared Error (TRMSE) with optional negative handling.
    
    Args:
        obs: Array of observed values
        sim: 2D array of simulated values [..., [year, month, day, flow]]
        handle_negatives: If True, uses abs() on both simulated and observed flows (original behavior)
                         If False, uses raw values (may detect negative flow issues)
    
    Returns:
        TRMSE value (float)
    """
    # Handle negatives based on parameter
    sim_flows = np.abs(sim) if handle_negatives else sim
    obs_flows = np.abs(obs) if handle_negatives else obs
    
    # Power transform
    z_sim = (np.power(sim_flows + 1, 0.3) - 1) / 0.3
    z_obs = (np.power(obs_flows + 1, 0.3) - 1) / 0.3
    
    return np.sqrt(np.mean((z_sim - z_obs) ** 2))

def calculate_msde(obs, sim, handle_negatives=True):
    """
    Calculate Mean Squared Derivative Error (MSDE) with optional negative handling.
    
    Args:
        obs: Array of observed values
        sim: 2D array of simulated values [..., [year, month, day, flow]]
        handle_negatives: If True, uses abs() on simulated flows (original behavior)
                         If False, uses raw values
    
    Returns:
        MSDE value (float)
    """
    sim_flows = np.abs(sim) if handle_negatives else sim
    
    obs_diff = np.diff(obs)
    sim_diff = np.diff(sim_flows)
    
    # Handle edge cases
    if len(obs_diff) == 0:
        return np.nan
    
    return np.mean((obs_diff - sim_diff) ** 2)

def calculate_roce(obs, sim, handle_negatives=True):
    """
    Calculate Runoff Coefficient Error (ROCE) with optional negative handling.
    
    Args:
        obs: Array of observed values
        sim: 2D array of simulated values [..., [year, month, day, flow]]
        handle_negatives: If True, uses abs() on simulated flows (original behavior)
                         If False, uses raw values
    
    Returns:
        ROCE value (float)
    """
    sim_flows = np.abs(sim) if handle_negatives else sim
    sum_sim = np.sum(sim_flows)
    sum_obs = np.sum(obs)
    
    # Handle division by zero
    if sum_obs == 0:
        return np.inf if sum_sim != 0 else np.nan
    
    return np.abs(sum_sim - sum_obs) / sum_obs
# %% ====================== RUNOFF MODEL FUNCTIONS ======================
def update_soil_parameters(variables,zone_id):
    """Update soil parameters for specific zone (zone_id) based on current variables"""
    global zones, dec_var,pop_ind,generation, vic_count
    
    # ===== 1. Load zone & soil files =====
    zone_lines = pd.read_csv(ZONE_FILE, delimiter='\t', header=None)
    zone_lines.columns = [f'Col{i+1}' for i in range(zone_lines.shape[1])]
    
    soil_lines = pd.read_csv(SOIL_FILE, delimiter='\t', header=None)
    soil_lines.columns = [f'Col{i+1}' for i in range(soil_lines.shape[1])]
    
    # ===== 2. Convert decision variables to array =====
    vars_array = np.array(variables, dtype=np.float64)
    var_lines = pd.DataFrame({'parameters': vars_array})
    
    # ===== 3. Update CURRENT ZONE parameters =====
    selected_values = var_lines.values.flatten()
    selected_values = [float(value) for value in selected_values]
    
    # Find coordinates for current zone
    filtered_rows = zone_lines[zone_lines['Col5'] == zone_id]
    filtered_coords = filtered_rows[['Col3', 'Col4']].drop_duplicates()
    soil_coords = soil_lines[['Col3', 'Col4']].reset_index()
    merged = pd.merge(soil_coords, filtered_coords, on=['Col3', 'Col4'], how='inner')
    matched_indices = merged['index']
    
    # Activate only grid cells in the current zone (set Col1 = 1); deactivate all others
    soil_lines['Col1'] = 0  # deactivate all
    soil_lines.loc[matched_indices, 'Col1'] = 1  # activate matching zone cells

    # Update soil parameters for current zone
    columns_to_update = ['Col5', 'Col6', 'Col7', 'Col8', 'Col9', 'Col19', 'Col20']
    for col, value in zip(columns_to_update, selected_values[:7]):
        soil_lines.loc[matched_indices, col] = round(float(value), 6)
    
    # ===== 4. Restore PREVIOUSLY CALIBRATED zones =====
    if zone_id > 1:
        for zone_id_prev in range(1,zone_id):
            # Load optimal parameters for previous zone
            param_file = os.path.join(TRACKING_DIRS['performance'], f'optimal_parameters_zone{zone_id_prev}.txt')
            optimal_param_prev=np.loadtxt(param_file)
            
            # Find coordinates for previous zone
            filtered_rows = zone_lines[zone_lines['Col5'] == zone_id_prev]
            filtered_coords = filtered_rows[['Col3', 'Col4']].drop_duplicates()
            soil_coords = soil_lines[['Col3', 'Col4']].reset_index()
            merged = pd.merge(soil_coords, filtered_coords, on=['Col3', 'Col4'], how='inner')
            matched_indices = merged['index']
            # Restore optimal parameters
            for col, value in zip(columns_to_update, optimal_param_prev[:7]):
                soil_lines.loc[matched_indices, col] = round(float(value), 6)
               
    
    # ===== 5. Save updated parameters =====
    soil_lines.to_csv(
        os.path.join(TRACKING_DIRS['soil'], 
        f'soilparam_zone{zone_id}_gen{generation}_pop{pop_ind}_round{vic_count}.txt'), 
        sep='\t', index=False, header=False
    )
    soil_lines.to_csv(SOIL_FILE, sep='\t', index=False, header=False)
    
    

def run_runoff_model():
    """Run the VIC rainfall-runoff model"""
    os.chdir(RUNOFF_SUBDIRS['output'])         
    os.chdir(RUNOFF_SUBDIRS['model'])
    os.system('./vicNl -g ../parameters/globalparam.txt')

# %% ====================== ROUTING MODEL FUNCTIONS ======================   
def update_routing_parameters(variables,zone_id):
    """Update routing parameters for specific zone (zone_id) based on current variables"""
    
    global generation, vic_count, dec_var, pop_ind
    global zones
    
    # ===== 1. Load zone & soil files =====
    zone_lines = pd.read_csv(ZONE_FILE, delimiter='\t', header=None)
    zone_lines.columns = [f'Col{i+1}' for i in range(zone_lines.shape[1])]
    
    # ===== 2. Convert decision variables to array =====
    vars_array = np.array(variables, dtype=np.float64)
    var_lines = pd.DataFrame({'parameters': vars_array}) 
    
    # ===== 3. Update CURRENT ZONE parameters =====
    selected_values = var_lines.values.flatten()
    selected_values = [float(value) for value in selected_values]
    # Extract velocity and diffusivity for current zone
    current_velocity = selected_values[-2]
    current_diffusivity = selected_values[-1]
    
    # Initialize maps with the current zone
    zone_velocity_map = {zone_id: current_velocity}
    zone_diffusivity_map = {zone_id: current_diffusivity}
    
    # ===== 4. Restore PREVIOUSLY CALIBRATED zones =====
    if zone_id > 1:
        for zone_id_prev in range(1, zone_id):
            # Load optimal parameters for previous zone
            param_file = os.path.join(TRACKING_DIRS['performance'], f'optimal_parameters_zone{zone_id_prev}.txt')
            optimal_param_prev = np.loadtxt(param_file)
            
            # Add previous zone's parameters to the maps
            zone_velocity_map[zone_id_prev] = optimal_param_prev[-2]
            zone_diffusivity_map[zone_id_prev] = optimal_param_prev[-1]
    
    # ===== 5. Apply all parameters to zone file =====
    zone_lines['velocity'] = zone_lines['Col5'].map(zone_velocity_map)
    zone_lines['diffusivity'] = zone_lines['Col5'].map(zone_diffusivity_map)
    
    # ===== 6. Save updated parameters =====     
    # Write routing parameters
    if dec_var[8] > 0 or dec_var[9] > 0:
        vel_ascii_path = os.path.join(
            ROUTING_SUBDIRS['parameters'], 
            f'velocity_zone{zone_id}_gen{generation}_pop{pop_ind}_round{vic_count}.txt')
        dif_ascii_path = os.path.join(
            ROUTING_SUBDIRS['parameters'], 
            f'diffusivity_zone{zone_id}_gen{generation}_pop{pop_ind}_round{vic_count}.txt')
        
        write_zone_param_ascii(zone_lines, ASCII_TEMP, vel_ascii_path, 'velocity')
        write_zone_param_ascii(zone_lines, ASCII_TEMP, dif_ascii_path, 'diffusivity')
        
    # ===== 7. Update config file =====
    with open(ROUTING_CONFIG, 'r') as f:
        base_config = f.read().splitlines()

    # Common parameter updates
    updated_lines = base_config.copy()
    if dec_var[8] > 0:  # Velocity
        updated_lines[4] = ".true."
        updated_lines[5] = f"../parameters/velocity_zone{zone_id}_gen{generation}_pop{pop_ind}_round{vic_count}.txt"
    if dec_var[9] > 0:  # Diffusivity
        updated_lines[7] = ".true."
        updated_lines[8] = f"../parameters/diffusivity_zone{zone_id}_gen{generation}_pop{pop_ind}_round{vic_count}.txt"

    # Single-station configurations
    station_lines = updated_lines.copy()
    station_lines[16] = f"../parameters/stations{zone_id}.txt"
    
    config_path = os.path.join(ROUTING_SUBDIRS['parameters'], f'config_station{zone_id}.txt')
    with open(config_path, 'w') as f:
        f.write("\n".join(station_lines[:42]))


def run_routing_model(zone_id):
    """
    Run routing model
    """
    os.chdir(ROUTING_SUBDIRS['model'])
    
    # Single-station mode - sequential processing
    config_path = f'../parameters/config_station{zone_id}.txt'
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Missing config: {config_path}")
    os.system(f'./rout {config_path}')

# %% ====================== OUTPUTS PROCESSING ======================  
def calculate_standardized_metrics(gauge_data, vic_data, sim_days):
    """Calculate all objective functions"""
    norm_trmse = 0
    norm_msde = 0
    mean_gauge_trmse = sum(gauge_data) / len(gauge_data)
    
    for i in range(len(gauge_data)):
        norm_msde += pow((gauge_data[i]-gauge_data[i-1])*2, 2)
        zobs = (pow((abs(gauge_data[i])+1), 0.3)-1)/0.3   # zobs is the power transformation
        norm_trmse += pow(zobs-mean_gauge_trmse, 2)
    
    nash = calculate_nse(gauge_data, vic_data)   
    trmse = calculate_trmse(gauge_data, vic_data)
    msde = calculate_msde(gauge_data, vic_data)
    roce = calculate_roce(gauge_data, vic_data)
    
    stand_nash = 1 if nash <= 0 else 1 - nash
    stand_trmse = min(max(trmse/pow(norm_trmse/(sim_days), 0.5), 0), 1)
    stand_msde = min(max(msde / norm_msde * (sim_days), 0), 1)
    stand_roce = min(max(roce, 0), 1)
    
    return {
        'nash': nash,
        'trmse': trmse,
        'msde': msde,
        'roce': roce,
        'norm_trmse': norm_trmse,
        'norm_msde': norm_msde,
        'stand_nash': stand_nash,
        'stand_trmse': stand_trmse,
        'stand_msde': stand_msde,
        'stand_roce': stand_roce
    }
  
def process_simulation_data(station,start_year_active, end_year,zone_id):
    """
    Aligns observed and simulated discharge data for a station.
    
    Args:
        station (str): Station name (without extension)
        start_year_active (int): First year after spinup
        end_year (int): Last simulation year
    
    Returns:
        tuple: (gauge_data, vic_data) - aligned arrays of observed and simulated flows
    """
    global generation, vic_count, pop_ind
    
    # 1. Source and destination paths
    station_name = f'{station}_zone{zone_id}_gen{generation}_pop{pop_ind}_round{vic_count}.day'
    src_path = os.path.join(ROUTING_SUBDIRS['output'], f'{station}.day')
    dest_path1 = os.path.join(ROUTING_SUBDIRS['output'], station_name)
    dest_path2 = os.path.join(TRACKING_DIRS['station'], station_name)
    # Perform the file copies
    shutil.copy(src_path, dest_path1)
    shutil.copy(src_path, dest_path2)
    
    # 2. Read observed discharge
    obs_file = os.path.join(ROUTING_SUBDIRS['obs_discharge'], f'{station}.csv')
    obs_df = pd.read_csv(obs_file)
    
    # Filter for active period
    obs_df = obs_df[(obs_df['Year'] >= start_year_active) & 
                    (obs_df['Year'] <= end_year)]
    gauge_data = obs_df['Flow'].values  
    
    # 3. Read simulated discharge 
    sim_file = os.path.join(ROUTING_SUBDIRS['output'], station_name)
    sim_data = []
    
    with open(sim_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 4:  # year, month, day, flow
                sim_data.append({
                    'year': int(parts[0]),
                    'month': int(parts[1]), 
                    'day': int(parts[2]),
                    'flow': float(parts[3])
                })
    
    # 4. Align data by date
    vic_data = []
    for _, obs_row in obs_df.iterrows():
        match = next((x for x in sim_data 
                     if x['year'] == obs_row['Year'] and
                     x['month'] == obs_row['Month'] and
                     x['day'] == obs_row['Day']), None)
        
        vic_data.append(match['flow'] if match else 0.0)
    
    return gauge_data, np.array(vic_data)
  
def process_calibration_objectives(site_names,zone_id):
    """Process results from one station (one zone) and return:
    - calibration_results: averaged standardized metrics for optimization
    - station_metrics: list of dicts with raw station metrics for tracking
    """
    global generation, vic_count, sim_period, hotstart, obj_fn,pop_ind
    
    
    # Initialize storage
    station_metrics = []
    
    # Get simulation years
    with open(ROUTING_CONFIG, 'r') as f:
        lines = f.read().splitlines()
    year_line = lines[30].strip()
    start_year = int(year_line.split()[0])
    end_year = int(year_line.split()[2])
    
    # Calculate active start year
    start_date = date(start_year, 1, 1)
    spinup_date = start_date + timedelta(days=hotstart)
    start_year_active = start_year + (spinup_date.year - start_date.year)
    sim_days= sim_period-hotstart
    
    
    station = site_names[zone_id-1][:-4]
    station_name = f'{station}_zone{zone_id}_gen{generation}_pop{pop_ind}_round{vic_count}.day'
    
    # Copy results
    shutil.copy(
        os.path.join(ROUTING_SUBDIRS['output'], f'{station}.day'),
        os.path.join(TRACKING_DIRS['station'], station_name)
    )
    
    # Process data
    gauge_data, vic_data = process_simulation_data(station, start_year_active, end_year,zone_id)
    results = calculate_standardized_metrics(gauge_data, vic_data,sim_days)
    
    # Store station metrics
    station_metrics.append({
        'station': station,
        'nash': results['nash'],
        'trmse': results['trmse'],
        'msde': results['msde'],
        'roce': results['roce'],
        'stand_nash': results['stand_nash'],
        'stand_trmse': results['stand_trmse'],
        'stand_msde': results['stand_msde'],
        'stand_roce': results['stand_roce']
    })
    
    
    
    # Prepare results
    calibration_results = []
    if obj_fn[0] > 0:
        calibration_results.append(results['stand_nash'])
        print(f"NSE = {results['nash']}")
    if obj_fn[1] > 0:
        calibration_results.append(results['stand_trmse'])
        print(f"TRMSE = {results['trmse']}")
    if len(obj_fn) > 2 and obj_fn[2] > 0:
        calibration_results.append(results['stand_msde'])
        print(f"MSDE = {results['msde']}")
    if len(obj_fn) > 3 and obj_fn[3] > 0:
        calibration_results.append(results['stand_roce'])
        print(f"ROCE = {results['roce']}")
    
    return calibration_results, station_metrics
# %% ====================== TRACKING PERFORMANCE METRICS ======================
def track_performance(variables, zone_id,station_metrics):
    """Track performance metrics and parameters"""
    global generation, vic_count,pop_ind,var_no
    
    # Create station-wise metrics DataFrame
    df_metrics_stations = pd.DataFrame(station_metrics)
    df_metrics_stations['generation'] = generation
    df_metrics_stations['population'] = pop_ind
    df_metrics_stations['vic_count'] = vic_count
    
    # Prepare parameters data
    vars_array = np.array(variables).reshape(var_no, 1, order='F')
    df_params_stations = pd.DataFrame(vars_array, columns=[f'station_{zone_id}'])
    df_params_stations['generation'] = generation
    df_params_stations['population'] = pop_ind
    df_params_stations['vic_count'] = vic_count
    
    # Write to files
    metrics_file = os.path.join(TRACKING_DIRS['performance'], 
                              f'track_metrics_zone{zone_id}_gen{generation}_pop{pop_ind}_round{vic_count}.txt')
    params_file = os.path.join(TRACKING_DIRS['performance'], 
                             f'track_parameters_zone{zone_id}_gen{generation}_pop{pop_ind}_round{vic_count}.txt')
    
    df_metrics_stations.to_csv(metrics_file, sep='\t', mode='a', index=False, header=not os.path.exists(metrics_file))
    df_params_stations.to_csv(params_file, sep='\t', mode='a', index=False, header=not os.path.exists(params_file))


# %% ====================== MAIN OPTIMIZATION FUNCTION or MAIN CALIBRATION MODEL======================
def viccall(variables):
    """Main calibration function (calling vic runoff/routing models) called by the optimization algorithm"""
    global vic_count, generation, pop_ind, multiprocessing_flag,zone_id
    global sim_period, hotstart, obj_fn, var_no, zones, dec_var,pop_size    
   
    # Calculate population index FIRST
    pop_ind = get_pop_ind(vic_count, generation, pop_size)
    
    if multiprocessing_flag == True:
        rank = multiprocessing.current_process()._identity[0]
        print(f"Thread no: {rank}, Round: {vic_count}") 
    else:
        print(f"Generation No: {generation}, Population Individual: {pop_ind}, Round: {vic_count}")
    
    print(f"Decision Variables [Soil/Routing Parameters]: {variables}")
    
    # List files and sort numerically based on the number in the filename (from upstream to downstream)
    site_names = sorted(os.listdir(ROUTING_SUBDIRS['obs_discharge']), key=lambda x: int(x[4:-4]))
    # site_name=os.path.splitext(site_names[zone_id-1])[0]
    
    # 1. Update all parameters
    update_soil_parameters(variables,zone_id)
    update_routing_parameters(variables,zone_id)
    
    # 2. Run models
    run_runoff_model()
    run_routing_model(zone_id)  # run for the station corresponding to the zone id
    
    # 3. Process and track results
    calibration_results,station_metrics = process_calibration_objectives(site_names,zone_id)
    track_performance(variables, zone_id, station_metrics)
    
    vic_count += 1
    return calibration_results

# %% ====================== CALIBRATION SETTINGS & RUN OPTIMIZATION ====================== 
def load_calibration_settings(filepath):
    """Load calibration settings from config file
    
    Returns:
        tuple: (sim_period, hotstart, iteration, pop_size, cores, 
               stations, zones, dec_var, obj_fn)
    """
    with open(filepath, 'r') as f:
        lines = f.read().split('\n')
    
    # Basic parameters
    sim_period = int(lines[2].split('\t')[0])
    hotstart = int(lines[3].split('\t')[0])
    iteration = int(lines[4].split('\t')[0])
    pop_size = int(lines[5].split('\t')[0])
    cores = int(lines[6].split('\t')[0])
    stations = int(lines[24].split('\t')[0])
    zones = int(lines[26].split('\t')[0])
    
    # Active parameters
    dec_var = [int(lines[i+8].split('\t')[0]) for i in range(10)]
    
    # Active objectives
    obj_fn = [int(lines[i+19].split('\t')[0]) for i in range(4)]
    
    return (sim_period, hotstart, iteration, pop_size, cores, 
            stations, zones, dec_var, obj_fn)

def setup_decision_variables(dec_var):
    """Initialize optimization parameter ranges"""
    param_ranges = {
        0: (0.000001, 0.9),      # b
        1: (0.000001, 0.999999),  # Ds
        2: (0.000001, 30),           # Dmax
        3: (0.000001, 0.999999),  # Ws
        4: (1, 3),            # c
        5: (0.05, 0.25),      # d1
        6: (0.3, 1.5),        # d2
        7: (1.5, 8.0),        # d3
        8: (0.5, 5),            # v
        9: (200, 4000)        # df
    }
    
    sel_var = []
    var_no = 0
    
    
    for i, active in enumerate(dec_var):
        if active == 1:
            var_no += 1
            low, high = param_ranges[i]
            sel_var.append(Real(low, high))

    return sel_var, var_no

def setup_objectives(obj_fn):
    """Configure optimization objectives"""
    min_obj = []
    max_obj = []
    eps = []
    obj_no = sum(obj_fn)
    
    for active in obj_fn:
        if active == 1:
            min_obj.append(0)
            max_obj.append(1)
            eps.append(0.001)
    
    return min_obj, max_obj, eps, obj_no 

def run_optimization(problem, eps, pop_size, max_iterations, cores):
    """
    Complete optimization workflow with execution and result saving in one function.
    
    Args:
        problem: Platypus Problem object
        eps: List of epsilon values for EpsNSGAII
        pop_size: Population size
        max_iterations: Maximum number of evaluations
        cores: Number of CPU cores for parallel processing
        
    Returns:
        tuple: (algorithm, duration) 
            algorithm: The optimization algorithm object
            duration: Total runtime as timedelta
    """
    global generation  # Declare we're using the global variable
    # Start timing
    start_time = datetime.datetime.now()
    
    # Initialize algorithm
    if multiprocessing_flag == True:
        with ProcessPoolEvaluator(cores) as evaluator:
            algorithm = EpsNSGAII(problem, eps, population_size=pop_size, evaluator=evaluator)
            # Run optimization loop
            while algorithm.nfe < max_iterations:
                print(f"Generation = {generation}")
                algorithm.step()
                generation += 1
    else:
        algorithm = EpsNSGAII(problem, eps, population_size=pop_size)
        # Run optimization loop
        while algorithm.nfe < max_iterations:
            print(f"Generation = {generation}")
            algorithm.step()
            generation += 1
    
    
    # Calculate total duration
    duration = datetime.datetime.now() - start_time
    
    return algorithm, duration


def save_nondominated_results(solutions,zone_id):
    """Save nondominated solutions to performance_tracking folder with metadata"""

    # Save nondominated results in text format
    params= [sol.variables[:] for sol in solutions]
    objectives= [sol.objectives[:] for sol in solutions]
    np.savetxt(
        os.path.join(TRACKING_DIRS['performance'], 'nondominated_parameters_zone%d.txt'%(zone_id)),
        params,
        fmt="%s"
    )
    np.savetxt(
        os.path.join(TRACKING_DIRS['performance'], 'nondominated_objectives_zone%d.txt'%(zone_id)),
        objectives,
        fmt="%s"
    )
    

def save_utopia_point_solution(solutions,zone_id):
    """
    Selects the optimal solution using Utopia point criterion.
    
    Parameters:
    -----------
    solutions : list
        List of solution objects from optimization (must have .objectives attribute)
        
    Returns:
    --------
    Solution
        The solution closest to the Utopia point (ideal objective vector)
        
    """
    # # 1. Get non-dominated front (already standardized 0-1)
    pareto_front = nondominated(solutions)
    
    if not pareto_front:
        raise ValueError("Empty Pareto front - no solutions to select from")
    
    # 2. Utopia point is [0, 0, ..., 0] for minimized [0,1] objectives
    objectives = np.array([sol.objectives for sol in pareto_front])
    distances = np.linalg.norm(objectives, axis=1)  # Euclidean distance to origin
    optimal_idx = np.argmin(distances)
    
    optimal_solution = pareto_front[optimal_idx]
    optimal_params = optimal_solution.variables
    optimal_objs = optimal_solution.objectives
    
    # 3. Save with zone-aware naming
    np.savetxt(
        os.path.join(TRACKING_DIRS['performance'], 'optimal_objectives_zone%d.txt'%(zone_id)),
        optimal_objs,
        fmt='%.6f'
    )
    
    np.savetxt(
        os.path.join(TRACKING_DIRS['performance'], 'optimal_parameters_zone%d.txt'%(zone_id)),
        optimal_params,
        fmt='%.6f'
    )
    
    return optimal_objs, optimal_params
# %% ====================== MAIN EXECUTION ======================
if __name__ == '__main__':   
    
    # 1. Initialize directories
    create_directories()
    validate_environment()
    
    # 2. Load calibration settings 
    (sim_period, hotstart, iteration, pop_size, cores, 
     stations, zones, dec_var, obj_fn) = load_calibration_settings(CALIB_SETTINGS)
    
    
    for zone_id in range(1,zones+1):
        print(f"\n{'='*50}")
        print(f"Starting optimization for Zone {zone_id}")
        print(f"{'='*50}")
        
        # 3. Set global variables explicitly
        vic_count = 1
        generation = 1
        pop_ind=1
        multiprocessing_flag = False
        #calibration_mode='zone'  # or 'basin'  [mean over zone or mean over basin]
        
        
        # 4. Setup optimization problem
        sel_var, var_no = setup_decision_variables(dec_var)
        min_obj, max_obj, eps, obj_no = setup_objectives(obj_fn)
        
        # 5. Initialize optimization problem
        problem = Problem(var_no, obj_no)
        problem.types[:] = sel_var
        problem.function = viccall
        hyp = Hypervolume(minimum=min_obj, maximum=max_obj)
        
        # 6. Run optimization
        print(f"Starting optimization. Population size: {pop_size}, Max evaluations: {iteration}")
        print(f"Active objectives: {obj_fn}")
        print(f"Decision variables: {dec_var}")
        
        algorithm, duration = run_optimization(
            problem, eps, pop_size, iteration, cores
        )
        
        # 7. Save final Pareto front
        nondominated_solutions = nondominated(algorithm.result)
        save_nondominated_results(nondominated_solutions,zone_id)
        print("\nOptimization complete!")
        print(f"Total time: {duration}")
        print(f"Nondominated solutions saved to: {TRACKING_DIRS['performance']}")
        
        
        # 8. Find & Save the Selected Optimal Solution using the Utopia Point
        save_utopia_point_solution(algorithm.result,zone_id)
        
        
        
        
    
   
   
    
    
    
   



## python /home/fs01/yl3984/VICRes/routing/parameters/mk_vicRes_input.py basin_name mode
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
import geopandas as gpd
from shapely.geometry import Point, Polygon
import sys
import os

case = sys.argv[1]      # basin5
sim_mode= sys.argv[2]   # water  # ener # cal # val
worker_id = sys.argv[3] # None
ystart  = 2010
yend    = 2020
if 'cal' or 'val' in sim_mode:
    inpath  = f'/home/fs01/yl3984/VIC/samples/benchmark/{worker_id}/'
    if 'cal' in sim_mode:
        outpath = '/home/fs01/yl3984/VICRes/routing/input/'+case+'/'
    else:
        gen_id  = int(worker_id.split('_')[0].replace('gen', ''))
        pop_id  = int(worker_id.split('_')[1].replace('pop', '')) 
        outpath = f'/home/fs01/yl3984/VICRes/routing/input/{case}/gen{gen_id}_pop{pop_id}/'
        if os.path.exists(outpath):
            sys.exit(0)
        os.makedirs(outpath, exist_ok=True)
else:
    inpath  = '/home/fs01/yl3984/VIC/samples/benchmark/10y_'+sim_mode+'/' 
    outpath = '/home/fs01/yl3984/VICRes/routing/input/'+case+'_ori/'
os.makedirs(outpath,exist_ok=True)
path_vic = '/home/fs01/yl3984/VIC/'

def cal_days():
    date1 = datetime.strptime(str(ystart)+"0101", "%Y%m%d")
    date2 = datetime.strptime(str(yend)+"1231", "%Y%m%d")
    delta = date2 - date1
    return int(delta.days+1)
period_days = cal_days()
print(period_days)

def read_nc(fname):
    nf = nc.Dataset(fname,'r')
    varname = nf.variables.keys()
    varname = list(varname)
    lat = np.array(nf.variables[varname[3]][:])
    lon = np.array(nf.variables[varname[2]][:])
    vic_trans_vege = np.array(nf.variables[varname[6]][:])
    runoff = np.array(nf.variables[varname[7]][:])
    baseflow = np.array(nf.variables[varname[8]][:])
    AIRTEMP   = np.array(nf.variables[varname[11]][:])
    RHD       = np.array(nf.variables[varname[17]][:])
    WINDS     = np.array(nf.variables[varname[18]][:])
    vic_eva_vege   = np.array(nf.variables[varname[19]][:])
    vic_eva_soil   = np.array(nf.variables[varname[20]][:])
    runoff = np.array(np.where(runoff<-10**3,np.nan,runoff))
    baseflow = np.array(np.where(baseflow<-10**3,np.nan,baseflow))
    return lat,lon,runoff,baseflow,vic_eva_vege,vic_trans_vege,vic_eva_soil,RHD,AIRTEMP,WINDS

lat_nc,lon_nc,_,_,_,_,_,_,_,_ = read_nc(inpath+'fluxes.'+str(ystart)+'-01-01.nc')
def read_VIC_output():
    output_all = np.full((period_days,len(lat_nc),len(lon_nc),8),np.nan)  # store runoff and baseflow
    day_index  = 0
    for yr in range(ystart,yend+1):
        _,_,runoff_each,baseflow_each,vic_eva_vege,vic_trans_vege,vic_eva_soil,RHD,AIRTEMP,WINDS = read_nc(inpath+'fluxes.'+str(yr)+'-01-01.nc')
        output_all[day_index:day_index+np.shape(runoff_each)[0],:,:,0] = runoff_each
        output_all[day_index:day_index+np.shape(runoff_each)[0],:,:,1] = baseflow_each
        output_all[day_index:day_index+np.shape(runoff_each)[0],:,:,2] = vic_eva_vege
        output_all[day_index:day_index+np.shape(runoff_each)[0],:,:,3] = vic_trans_vege
        output_all[day_index:day_index+np.shape(runoff_each)[0],:,:,4] = vic_eva_soil
        output_all[day_index:day_index+np.shape(runoff_each)[0],:,:,5] = RHD
        output_all[day_index:day_index+np.shape(runoff_each)[0],:,:,6] = AIRTEMP
        output_all[day_index:day_index+np.shape(runoff_each)[0],:,:,7] = WINDS
        day_index = day_index + np.shape(runoff_each)[0]
    return output_all
output_all = read_VIC_output()

def basin_region(case):
    region = np.array([[-124, -120, 38.5, 42.5],
                             [-120, -114.5, 38, 42],
                             [-114, -109, 41, 45],
                             [-113, -107, 30, 35],
                             [-110, -105, 37.5, 42],
                             [-107, -101, 32, 37.5],
                             [-103, -96.5, 37.5, 40],
                             [-100.54, -92, 42.5, 48],
                             [-99, -94, 30, 35],
                             [-89, -83, 38, 42.5],
                             [-86, -83, 29, 36],
                             [-80, -74.5, 40, 43],
                             [-74, -70, 41, 45.5]])
    return region[int(case[5:]),0],region[int(case[5:]),1],region[int(case[5:]),2],region[int(case[5:]),3]
lon1, lon2, lat1, lat2 = basin_region(case)
shpfile = path_vic + "/samples/basin/Basins_NSF_definitive.shp"

def seperate_connected_basin(lon1, lon2, lat1, lat2, shpfile):
    # Read shapefile
    gdf = gpd.read_file(shpfile)
    # Convert geometry to lon/lat
    if gdf.crs is None:
        gdf.set_crs("EPSG:5070", inplace=True)
    # Convert to lat/lon
    gdf = gdf.to_crs("EPSG:4326")
    # Select basins by bounding box
    gdf_select = gdf.cx[lon1:lon2, lat1:lat2]
    gdf_exploded = gdf_select.explode(index_parts=True, ignore_index=True)
    print(len(gdf_exploded))
    gdf_plot = gdf_exploded.iloc[[0]]  # iloc=0, basin5; iloc=1, basin6
    return gdf_plot
gdf_plot = seperate_connected_basin(lon1, lon2, lat1, lat2, shpfile)

def basin_grids(gdf_basin, lon_result, lat_result):
    polygon = gdf_basin.geometry.iloc[0]
    # Extract polygon bounds
    lon_poly, lat_poly = polygon.exterior.coords.xy
    lon_min, lon_max = np.min(lon_poly), np.max(lon_poly)
    lat_min, lat_max = np.min(lat_poly), np.max(lat_poly)
    # Create masks for 1D arrays
    mask_lon = (lon_result >= lon_min) & (lon_result <= lon_max)
    mask_lat = (lat_result >= lat_min) & (lat_result <= lat_max)
    # Combine to get 2D grid points inside bounding box
    mask_2d = np.outer(mask_lat, mask_lon)  # shape (lat, lon)
    # Indices of points inside bounding box
    row_inside, col_inside = np.where(mask_2d)
    # Extract coordinates
    lon_point = lon_result[col_inside]
    lat_point = lat_result[row_inside]
    return lon_point, lat_point, row_inside, col_inside
lon_basin_txt, lat_basin_txt, row_basin_txt, col_basin_txt = basin_grids(gdf_plot, lon_nc, lat_nc)

def write_txt():
    for file_ind in range(len(lon_basin_txt)):
        lon_txt = lon_basin_txt[file_ind]
        lat_txt = lat_basin_txt[file_ind]
        # Extract runoff and baseflow time series for this grid
        runoff_txt = output_all[:, row_basin_txt[file_ind], col_basin_txt[file_ind], 0]
        base_txt   = output_all[:, row_basin_txt[file_ind], col_basin_txt[file_ind], 1]
        vic_evege_txt   = output_all[:, row_basin_txt[file_ind], col_basin_txt[file_ind], 2]
        vic_tvege_txt   = output_all[:, row_basin_txt[file_ind], col_basin_txt[file_ind], 3]
        vic_esoil_txt   = output_all[:, row_basin_txt[file_ind], col_basin_txt[file_ind], 4]
        RHD_txt   = output_all[:, row_basin_txt[file_ind], col_basin_txt[file_ind], 5]*100
        AIRTEMP_txt = output_all[:, row_basin_txt[file_ind], col_basin_txt[file_ind], 6]
        WINDS_txt   = output_all[:, row_basin_txt[file_ind], col_basin_txt[file_ind], 7]
        # Output file name
        filename = f"fluxes_{lat_txt:.4f}_{lon_txt:.4f}"
        # Start date
        date1 = datetime.strptime(f"{ystart}0101", "%Y%m%d")
        # Write to file
        with open(outpath + filename, 'w') as f:
            for day_ind in range(period_days):
                date_c = date1 + timedelta(days=day_ind)
                yyyy = date_c.year
                mm = date_c.month
                dd = date_c.day
                hh = 0  # fixed to 0 hour
                # Convert to string with 4-decimal precision
                runoff_str = f"{runoff_txt[day_ind]:.4f}"
                base_str   = f"{base_txt[day_ind]:.4f}"
                vic_evege_str   = f"{vic_evege_txt[day_ind]:.4f}"
                vic_tvege_str   = f"{vic_tvege_txt[day_ind]:.4f}"
                vic_esoil_str   = f"{vic_esoil_txt[day_ind]:.4f}"
                RHD_str   = f"{RHD_txt[day_ind]:.4f}"
                AIRTEMP_str = f"{AIRTEMP_txt[day_ind]:.4f}"
                WINDS_str   = f"{WINDS_txt[day_ind]:.4f}"
                # Write line to file
                f.write(f"{yyyy}\t{mm}\t{dd}\t{hh}\t{runoff_str}\t{base_str}\t{vic_evege_str}\t{vic_tvege_str}\t{vic_esoil_str}\t{RHD_str}\t{AIRTEMP_str}\t{WINDS_str}\n")
            f.close()
write_txt()

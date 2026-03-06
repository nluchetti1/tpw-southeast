import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy import Scene
import metpy.calc as mpcalc
from metpy.units import units
from herbie import Herbie
import datetime
import xarray as xr
import numpy as np
import requests
import os
import imageio.v2 as imageio

# --- SETUP ---
now = datetime.datetime.utcnow().replace(minute=0, second=0, microsecond=0)
frames = []
num_steps = 8  

if not os.path.exists('frames'):
    os.makedirs('frames')

for i in range(num_steps, -1, -1):
    target_time = now - datetime.timedelta(hours=i)
    time_str = target_time.strftime("%Y%m%d_%H00")
    
    # ... [Satellite code remains same] ...

    # 2. FETCH RAP WITH VARIABLE INSPECTOR
    try:
        H = Herbie(target_time, model='rap', product='awp130pgrb', fxx=0, 
                   priority=['aws', 'nomads'], verbose=False)
        
        # Load the dataset
        ds_rap = H.xarray(":(SPFH|UGRD|VGRD):925 mb").metpy.parse_cf()
        
        # --- VARIABLE INSPECTOR PRINT BLOCK ---
        print(f"\n--- INSPECTING RAP {time_str}Z ---")
        for var in ds_rap.data_vars:
            long_name = ds_rap[var].attrs.get('long_name', 'No description')
            print(f"Variable Key: '{var}' | Description: {long_name}")
        print("-----------------------------------\n")
        
        # Logic to handle whatever name it finds
        u = ds_rap['u'].metpy.unit_array if 'u' in ds_rap.data_vars else ds_rap['ugrd'].metpy.unit_array
        v = ds_rap['v'].metpy.unit_array if 'v' in ds_rap.data_vars else ds_rap['vgrd'].metpy.unit_array
        
        # Check for 'q' or 'spfh' based on common GRIB2 mapping
        if 'q' in ds_rap.data_vars:
            q = ds_rap['q'].metpy.unit_array
        elif 'spfh' in ds_rap.data_vars:
            q = ds_rap['spfh'].metpy.unit_array
        else:
            raise KeyError("Could not find humidity variable in dataset")
            
        lons, lats = ds_rap.longitude, ds_rap.latitude
        dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
        adv = mpcalc.advection(q, u, v, dx=dx, dy=dy) * 1e6
        
    except Exception as e:
        print(f"  > RAP extraction failed for {time_str}: {e}")

    # ... [Plotting code remains same] ...

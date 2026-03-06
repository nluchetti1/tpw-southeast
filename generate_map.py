import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy import Scene
import metpy.calc as mpcalc
from metpy.units import units
from siphon.catalog import TDSCatalog
import datetime
import xarray as xr
import numpy as np
import requests
import os

# --- 1. DOWNLOAD THE LATEST MIMIC NETCDF ---
# We check the last 4 hours of the SSEC archive to find the latest valid file
now = datetime.datetime.utcnow()
local_file = "latest_mimic.nc"
found_file = False

for i in range(4):
    check_time = now - datetime.timedelta(hours=i)
    # File naming convention: YYYYMMDD_HH00.nc
    fname = check_time.strftime("%Y%m%d_%H00.nc")
    url = f"https://tropic.ssec.wisc.edu/archive/data/mtpw2/{check_time.year}/{fname}"
    
    try:
        # SSEC often requires a standard User-Agent to prevent 403/404 errors
        headers = {'User-Agent': 'Mozilla/5.0'}
        r = requests.get(url, headers=headers, timeout=15)
        if r.status_code == 200:
            with open(local_file, 'wb') as f:
                f.write(r.content)
            found_file = True
            print(f"Successfully downloaded: {fname}")
            break
    except Exception as e:
        print(f"Trying {fname} failed: {e}")

# --- 2. PROCESS WITH SATPY ---
if found_file:
    try:
        # Load the local file into Satpy
        scn = Scene(reader='mimic_TPW2_nc', filenames=[local_file])
        scn.load(['tpw'])
        
        # Convert to xarray and handle potential y-axis flipping
        tpw_ds = scn.to_xarray_dataset()
        tpw_data = tpw_ds['tpw']
        if tpw_data.y[0] < tpw_data.y[-1]:
            tpw_data = tpw_data.sortby('y', ascending=False)
            
        # Subset to the Southeast US
        tpw_data = tpw_data.sel(y=slice(42, 28), x=slice(-90, -70))
    except Exception as e:
        print(f"Satpy Processing Error: {e}")
        tpw_data = None
else:
    tpw_data = None

# --- 3. FETCH RAP 13KM (925MB) ---
try:
    cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml')
    ds_rap = cat.datasets[0].remote_access(use_xarray=True)
    
    # RAP uses 'isobaric' for pressure levels
    subset_rap = ds_rap.metpy.sel(isobaric=925 * units.hPa, 
                                 lat=slice(42, 28), lon=slice(360-90, 360-70))
    
    u = subset_rap['u-component_of_wind_isobaric'].metpy.unit_array
    v = subset_rap['v-component_of_wind_isobaric'].metpy.unit_array
    q = subset_rap['Specific_humidity_isobaric'].metpy.unit_array
    lats_r, lons_r = subset_rap.lat, subset_rap.lon
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lons_r, lats_r)
    adv = mpcalc.advection(q, u, v, dx=dx, dy=dy) * 1e6 # Scaled for contouring
except Exception as e:
    print(f"RAP/MetPy Error: {e}")
    adv = None

# --- 4. COMPOSITE PLOT ---
fig = plt.figure(figsize=(14, 10), facecolor='black')
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([-88, -74, 31, 40]) # NC/VA/SC Sector

if tpw_data is not None:
    im = ax.pcolormesh(tpw_data.x, tpw_data.y, tpw_data, cmap='viridis', alpha=0.8, shading='auto')
    cb = plt.colorbar(im, ax=ax, pad=0.02, aspect=30)
    cb.set_label('Total Precipitable Water (mm)', color='white')
    cb.ax.yaxis.set_tick_params(color='white', labelcolor='white')

if adv is not None:
    # Overlay RAP Moisture Advection (Positive = Red)
    clevs = np.arange(2, 22, 4)
    cs = ax.contour(lons_r, lats_r, adv, levels=clevs, colors='red', linewidths=1.5)
    ax.clabel(cs, inline=True, fontsize=10, fmt='%d', colors='white')
    
    # Overlay 925mb Wind Vectors
    ax.quiver(lons_r[::3], lats_r[::3], u[::3, ::3], v[::3, ::3], 
              color='white', scale=400, width=0.003, alpha=0.9)

ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='white', linewidth=1)
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), edgecolor='cyan')

plt.title(f"MIMIC-TPW2 + RAP 925mb Moisture Advection\nValid: {datetime.datetime.utcnow().strftime('%H:%MZ %d %b %Y')}", 
          color='white', fontsize=14, pad=20)

plt.savefig('output_map.png', facecolor='black', bbox_inches='tight', dpi=150)

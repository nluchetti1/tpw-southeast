import matplotlib
matplotlib.use('Agg')
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

# --- 1. SATELLITE DATA (MIMIC-TPW2) ---
now = datetime.datetime.utcnow()
local_file = "latest_mimic.nc"
found_sat = False

for i in range(4):
    check_time = now - datetime.timedelta(hours=i)
    fname = check_time.strftime("%Y%m%d_%H00.nc")
    url = f"https://tropic.ssec.wisc.edu/archive/data/mtpw2/{check_time.year}/{fname}"
    try:
        r = requests.get(url, headers={'User-Agent': 'Mozilla/5.0'}, timeout=15)
        if r.status_code == 200:
            with open(local_file, 'wb') as f:
                f.write(r.content)
            found_sat = True
            break
    except: continue

tpw_data = None
if found_sat:
    try:
        scn = Scene(reader='mimic_TPW2_nc', filenames=[local_file])
        scn.load(['tpw'])
        tpw_ds = scn.to_xarray_dataset()
        tpw_data = tpw_ds['tpw']
        if tpw_data.y[0] < tpw_data.y[-1]:
            tpw_data = tpw_data.sortby('y', ascending=False)
        tpw_data = tpw_data.sel(y=slice(42, 28), x=slice(-90, -70))
    except Exception as e: print(f"Satpy Error: {e}")

# --- 2. MODEL DATA (RAP 13km) ---
adv = None
try:
    cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml')
    ds_rap = cat.datasets[0].remote_access(use_xarray=True).metpy.parse_cf()
    
    # Selecting 925mb. We use MetPy selectors to handle the projected grid
    # 'isobaric' is the vertical coordinate in RAP
    subset_rap = ds_rap.metpy.sel(isobaric=925 * units.hPa, 
                                 lat=slice(42, 28), lon=slice(-90, -70))
    
    u = subset_rap['u-component_of_wind_isobaric']
    v = subset_rap['v-component_of_wind_isobaric']
    q = subset_rap['Specific_humidity_isobaric']
    
    # FORCE COORDINATE EXTRACTION: 
    # If standard 'lat' fails, we pull from the CF-parsed coordinates
    lons = subset_rap.metpy.longitude
    lats = subset_rap.metpy.latitude
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    adv = mpcalc.advection(q, u, v, dx=dx, dy=dy) * 1e6
except Exception as e: 
    print(f"RAP Error: {e}")

# --- 3. PLOT ---
fig = plt.figure(figsize=(12, 9), facecolor='black')
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([-88, -74, 31, 40]) # NC/VA/SC/GA Focus

# Background Satellite
if tpw_data is not None:
    sat_plot = ax.pcolormesh(tpw_data.x, tpw_data.y, tpw_data, cmap='viridis', alpha=0.8, shading='auto')
    plt.colorbar(sat_plot, ax=ax, orientation='vertical', pad=0.02, label='TPW (mm)')

# Foreground RAP Advection & Winds
if adv is not None:
    # Plot Advection Contours
    clevs = np.arange(2, 22, 4)
    cs = ax.contour(lons, lats, adv, levels=clevs, colors='red', linewidths=1.5, transform=ccrs.PlateCarree())
    ax.clabel(cs, inline=True, fontsize=10, fmt='%d', colors='white')
    
    # Plot Wind Quivers (Winds at 925mb)
    # Subsample [::4] to keep the map from getting cluttered
    ax.quiver(lons.values[::4, ::4], lats.values[::4, ::4], 
              u.values[::4, ::4], v.values[::4, ::4], 
              color='white', scale=400, transform=ccrs.PlateCarree())

# Geography
ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='white', linewidth=1)
ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='cyan')

plt.title(f"MIMIC-TPW2 Satellite + RAP 925mb Advection\nValid: {now.strftime('%H:%MZ %d %b %Y')}", 
          color='white', fontsize=14, pad=20)

plt.savefig('output_map.png', facecolor='black', bbox_inches='tight', dpi=150)

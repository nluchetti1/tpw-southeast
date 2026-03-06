import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy import Scene
import metpy.calc as mpcalc
from metpy.units import units
from herbie import Herbie # Use Herbie like your other script
import datetime
import xarray as xr
import numpy as np
import requests

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

# --- 2. THE HERBIE FIX FOR RAP 925MB ---
adv = None
lons, lats = None, None
u, v = None, None

try:
    # Look for the freshest RAP run (matching your other script's logic)
    for h_back in range(0, 4):
        check_time = (now - datetime.timedelta(hours=h_back)).replace(minute=0, second=0, microsecond=0)
        try:
            # Targeting the RAP model at 925mb for winds and moisture
            H = Herbie(check_time, model='rap', product='awp130pgrb', verbose=False)
            # Fetch 925mb Winds and Specific Humidity (used for moisture advection)
            ds_rap = H.xarray(":(UGRD|VGRD|SPFH):925 mb").metpy.parse_cf()
            break
        except: continue

    # Extract components using the names Herbie standardizes
    u = ds_rap['u'].metpy.unit_array
    v = ds_rap['v'].metpy.unit_array
    q = ds_rap['spfh'].metpy.unit_array # Specific Humidity at 925mb
    
    lats = ds_rap.latitude
    lons = ds_rap.longitude
    
    # Calculate Advection
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    adv = mpcalc.advection(q, u, v, dx=dx, dy=dy) * 1e6
except Exception as e: 
    print(f"Herbie RAP Extraction Failed: {e}")

# --- 3. PLOT ---
fig = plt.figure(figsize=(12, 9), facecolor='black')
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([-88, -74, 31, 40]) 

if tpw_data is not None:
    ax.pcolormesh(tpw_data.x, tpw_data.y, tpw_data, cmap='viridis', alpha=0.8, shading='auto')

if adv is not None:
    cs = ax.contour(lons, lats, adv, levels=np.arange(2, 22, 4), 
                    colors='red', linewidths=1.5, transform=ccrs.PlateCarree())
    ax.clabel(cs, inline=True, fontsize=10, fmt='%d', colors='white')
    ax.quiver(lons.values[::4, ::4], lats.values[::4, ::4], 
              u.values[::4, ::4], v.values[::4, ::4], 
              color='white', scale=400, transform=ccrs.PlateCarree())

ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='white', linewidth=0.8)
ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='cyan')

plt.title(f"MIMIC-TPW2 + RAP 925mb Moisture Advection\nValid: {now.strftime('%H:%MZ %d %b %Y')}", 
          color='white', fontsize=14, pad=15)

plt.savefig('output_map.png', facecolor='black', bbox_inches='tight', dpi=150)

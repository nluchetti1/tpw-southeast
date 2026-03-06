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
        # Conversion: mm to inches
        tpw_data = tpw_ds['tpw'] / 25.4 
        
        if tpw_data.y[0] < tpw_data.y[-1]:
            tpw_data = tpw_data.sortby('y', ascending=False)
        tpw_data = tpw_data.sel(y=slice(42, 28), x=slice(-90, -70))
    except Exception as e: print(f"Satpy Error: {e}")

# --- 2. MODEL DATA (RAP 925mb via Herbie) ---
adv = None
try:
    H = Herbie(now.replace(minute=0, second=0, microsecond=0), 
               model='rap', product='awp130pgrb', fxx=0, verbose=False)
    
    ds_rap = H.xarray(":(SPFH|UGRD|VGRD):925 mb").metpy.parse_cf()
    
    u = ds_rap['u'].metpy.unit_array
    v = ds_rap['v'].metpy.unit_array
    
    # Handle variable naming (q or spfh)
    q_key = 'q' if 'q' in ds_rap.data_vars else 'spfh'
    q = ds_rap[q_key].metpy.unit_array
    
    lons, lats = ds_rap.longitude, ds_rap.latitude
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    adv = mpcalc.advection(q, u, v, dx=dx, dy=dy) * 1e6
except Exception as e:
    print(f"Herbie RAP Error: {e}")

# --- 3. PLOT ---
fig = plt.figure(figsize=(12, 9), facecolor='black')
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([-88, -74, 31, 40]) 

# Background Satellite (Now in Inches)
if tpw_data is not None:
    # Typical TPW range in inches is 0.5 to 2.5
    sat_plot = ax.pcolormesh(tpw_data.x, tpw_data.y, tpw_data, 
                             cmap='viridis', alpha=0.7, shading='auto', vmin=0.5, vmax=2.5)
    cb = plt.colorbar(sat_plot, ax=ax, orientation='vertical', pad=0.02, aspect=30)
    cb.set_label('Precipitable Water (inches)', color='white')
    cb.ax.yaxis.set_tick_params(color='white', labelcolor='white')

# RAP Advection Overlay
if adv is not None:
    clevs = np.arange(2, 22, 4)
    cs = ax.contour(lons, lats, adv, levels=clevs, colors='red', linewidths=1.2, transform=ccrs.PlateCarree())
    ax.clabel(cs, inline=True, fontsize=9, fmt='%d', colors='white')
    ax.quiver(lons.values[::4, ::4], lats.values[::4, ::4], 
              u.values[::4, ::4], v.values[::4, ::4], 
              color='white', scale=400, transform=ccrs.PlateCarree(), alpha=0.8)

# Geography Features (with 10m Counties)
ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='white', linewidth=1.2)
ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='cyan', linewidth=1.0)
counties = cfeature.NaturalEarthFeature(category='cultural', name='admin_2_counties', 
                                        scale='10m', facecolor='none')
ax.add_feature(counties, edgecolor='gray', linewidth=0.4, alpha=0.5)

plt.title(f"MIMIC-TPW2 (in) + RAP 925mb Moisture Advection\nValid: {now.strftime('%H:%MZ %d %b %Y')}", 
          color='white', fontsize=14, pad=15)

plt.savefig('output_map.png', facecolor='black', bbox_inches='tight', dpi=150)

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
    
    tpw_data = None
    adv, lons, lats, u, v = None, None, None, None, None

    # 1. FETCH SATELLITE (Inches)
    fname = f"{time_str}.nc"
    url = f"https://tropic.ssec.wisc.edu/archive/data/mtpw2/{target_time.year}/{fname}"
    local_nc = f"sat_{i}.nc"
    
    try:
        r = requests.get(url, headers={'User-Agent': 'Mozilla/5.0'}, timeout=10)
        if r.status_code == 200:
            with open(local_nc, 'wb') as f:
                f.write(r.content)
            scn = Scene(reader='mimic_TPW2_nc', filenames=[local_nc])
            scn.load(['tpw'])
            ds_sat = scn.to_xarray_dataset()
            tpw_data = ds_sat['tpw'] / 25.4 # Convert mm to inches
            if tpw_data.y[0] < tpw_data.y[-1]:
                tpw_data = tpw_data.sortby('y', ascending=False)
            tpw_data = tpw_data.sel(y=slice(42, 28), x=slice(-90, -70))
    except: pass

    # 2. FETCH RAP (Using 'conus' for full variables)
    try:
        H = Herbie(target_time, model='rap', product='conus', fxx=0, 
                   priority=['aws', 'nomads'], verbose=False)
        
        # Pull 925mb level. 'conus' product includes humidity.
        ds_rap = H.xarray(":925 mb").metpy.parse_cf()
        
        u, v = ds_rap['u'].metpy.unit_array, ds_rap['v'].metpy.unit_array
        
        # Map humidity (standard name in RAP conus is 'q')
        q_key = 'q' if 'q' in ds_rap.data_vars else 'spfh'
        q = ds_rap[q_key].metpy.unit_array
            
        lons, lats = ds_rap.longitude, ds_rap.latitude
        dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
        adv = mpcalc.advection(q, u, v, dx=dx, dy=dy) * 1e6
    except Exception as e:
        print(f"  > RAP extraction failed for {time_str}: {e}")

    # 3. PLOT FRAME
    if tpw_data is not None or adv is not None:
        fig = plt.figure(figsize=(12, 9), facecolor='black')
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.set_extent([-88, -74, 31, 40])

        if tpw_data is not None:
            im = ax.pcolormesh(tpw_data.x, tpw_data.y, tpw_data, cmap='viridis', 
                               alpha=0.8, shading='auto', vmin=0.5, vmax=2.5)
            # Add colorbar for satellite data
            cb = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02, aspect=30)
            cb.set_label('Precipitable Water (inches)', color='white')
            cb.ax.yaxis.set_tick_params(color='white', labelcolor='white')

        if adv is not None:
            cs = ax.contour(lons, lats, adv, levels=np.arange(2, 22, 4), colors='red', linewidths=1.5)
            ax.clabel(cs, inline=True, fontsize=10, fmt='%d', colors='white')
            ax.quiver(lons.values[::4, ::4], lats.values[::4, ::4], u.values[::4, ::4], v.values[::4, ::4], 
                      color='white', scale=400, transform=ccrs.PlateCarree())

        # Features
        ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='white', linewidth=1.2)
        ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='cyan')
        counties = cfeature.NaturalEarthFeature('cultural', 'admin_2_counties', '10m', facecolor='none')
        ax.add_feature(counties, edgecolor='gray', linewidth=0.4, alpha=0.5)

        plt.title(f"MIMIC-TPW (in) + RAP 925mb Advection\nValid: {target_time.strftime('%H:%MZ %d %b %Y')}", 
                  color='white', fontsize=14)
        
        frame_path = f"frames/frame_{i:02d}.png"
        plt.savefig(frame_path, facecolor='black', bbox_inches='tight', dpi=120)
        frames.append(imageio.imread(frame_path))
        plt.close()
    
    if os.path.exists(local_nc):
        os.remove(local_nc)

# 4. SAVE GIF
if frames:
    imageio.mimsave('output_animation.gif', frames, fps=2, loop=0)
    imageio.imwrite('output_map.png', frames[-1])

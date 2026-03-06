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
num_steps = 10 

if not os.path.exists('frames'):
    os.makedirs('frames')

for i in range(num_steps, -1, -1):
    target_time = now - datetime.timedelta(hours=i)
    time_str = target_time.strftime("%Y%m%d_%H00")
    
    tpw_data = None
    transport, lons, lats, u, v = None, None, None, None, None

    # 1. FETCH SATELLITE (Inches)
    fname = f"{time_str}.nc"
    url = f"https://tropic.ssec.wisc.edu/archive/data/mtpw2/{target_time.year}/{fname}"
    local_nc = f"sat_{i}.nc"
    try:
        r_sat = requests.get(url, headers={'User-Agent': 'Mozilla/5.0'}, timeout=10)
        if r_sat.status_code == 200:
            with open(local_nc, 'wb') as f:
                f.write(r_sat.content)
            scn = Scene(reader='mimic_TPW2_nc', filenames=[local_nc])
            scn.load(['tpw'])
            ds_sat = scn.to_xarray_dataset()
            tpw_data = ds_sat['tpw'] / 25.4 
            if tpw_data.y[0] < tpw_data.y[-1]:
                tpw_data = tpw_data.sortby('y', ascending=False)
            tpw_data = tpw_data.sel(y=slice(42, 28), x=slice(-90, -70))
    except: pass

    # 2. FETCH RAP & CALCULATE MOISTURE TRANSPORT
    try:
        H = Herbie(target_time, model='rap', product='awp252pgrb', fxx=0, 
                   priority=['aws', 'nomads'], verbose=False)
        
        ds_rap = H.xarray(":925 mb").metpy.parse_cf()
        
        u_qty = ds_rap['u'].metpy.unit_array
        v_qty = ds_rap['v'].metpy.unit_array
        wind_speed = mpcalc.wind_speed(u_qty, v_qty)
        
        # Calculate Mixing Ratio (g/g) from T and RH
        rel_hum = ds_rap['r'].metpy.unit_array / 100.0
        temp = ds_rap['t'].metpy.unit_array
        pressure = 925 * units.hPa
        mixing_ratio = mpcalc.mixing_ratio_from_relative_humidity(pressure, temp, rel_hum)
        
        # Scale moisture transport by 100 as requested
        # We use .magnitude because contourf/quiver need raw arrays
        transport = (wind_speed.magnitude * mixing_ratio.magnitude) * 100
        u = u_qty.magnitude
        v = v_qty.magnitude
        lons = ds_rap.longitude.values
        lats = ds_rap.latitude.values
            
    except Exception as e:
        print(f"  > RAP process failed for {time_str}: {e}")

    # 3. PLOT FRAME
    if tpw_data is not None or transport is not None:
        fig = plt.figure(figsize=(12, 9), facecolor='black')
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.set_extent([-88, -74, 31, 40])

        if tpw_data is not None:
            # Satellite background in grayscale to let transport colors pop
            ax.pcolormesh(tpw_data.x, tpw_data.y, tpw_data, cmap='Greys_r', 
                          alpha=0.3, shading='auto', vmin=0.5, vmax=2.5)

        if transport is not None:
            # Pink shades for high transport values
            clevs = [5, 10, 15, 20, 24, 28, 32, 36, 40]
            cf = ax.contourf(lons, lats, transport, levels=clevs, cmap='RdPu', alpha=0.7)
            cb = plt.colorbar(cf, ax=ax, orientation='vertical', pad=0.02, aspect=30)
            cb.set_label('925mb Moisture Transport (m/s)', color='white')
            cb.ax.yaxis.set_tick_params(color='white', labelcolor='white')
            
            # Wind quivers - using .magnitude arrays for plotting
            ax.quiver(lons[::4, ::4], lats[::4, ::4], u[::4, ::4], v[::4, ::4], 
                      color='white', scale=400, transform=ccrs.PlateCarree(), alpha=0.6)

        ax.add_feature(cfeature.STATES.with_scale('10m'), edgecolor='white', linewidth=1.2)
        ax.add_feature(cfeature.COASTLINE.with_scale('10m'), edgecolor='cyan')
        counties = cfeature.NaturalEarthFeature('cultural', 'admin_2_counties', '10m', facecolor='none')
        ax.add_feature(counties, edgecolor='gray', linewidth=0.4, alpha=0.5)

        plt.title(f"MIMIC-TPW + RAP 925mb Moisture Transport\nValid: {target_time.strftime('%H:%MZ %d %b %Y')}", 
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

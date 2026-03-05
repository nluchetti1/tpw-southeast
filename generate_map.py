import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
from siphon.catalog import TDSCatalog
import datetime
import numpy as np

# 1. Fetch RAP 13km Data (925mb)
try:
    rap_cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml')
    ds_rap = rap_cat.datasets[0].remote_access(use_xarray=True)
    subset_rap = ds_rap.metpy.sel(vertical=925 * units.hPa, 
                                 lat=slice(42, 28), lon=slice(360-90, 360-70))
    
    u = subset_rap['u-component_of_wind_isobaric'].metpy.unit_array
    v = subset_rap['v-component_of_wind_isobaric'].metpy.unit_array
    q = subset_rap['Specific_humidity_isobaric'].metpy.unit_array
    lats_r, lons_r = subset_rap.lat, subset_rap.lon
    
    # Calculate Advection
    dx, dy = mpcalc.lat_lon_grid_deltas(lons_r, lats_r)
    adv = mpcalc.advection(q, u, v, dx=dx, dy=dy) * 1e6 
except Exception as e:
    print(f"RAP Access Error: {e}")

# 2. Fetch MIMIC-TPW2 Satellite Data (netCDF)
tpw_url = "https://tropic.ssec.wisc.edu/real-time/mtpw2/data/current.nc"
ds_tpw = xr.open_dataset(tpw_url)
tpw_data = ds_tpw['tpw'].sel(lat=slice(28, 42), lon=slice(-90, -70))

# 3. Create the Composite Plot
fig = plt.figure(figsize=(14, 9), facecolor='black')
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([-88, -74, 31, 40]) # Zoom on NC/SC/VA/GA

# Layer 1: Satellite TPW with Colorbar
im = ax.pcolormesh(tpw_data.lon, tpw_data.lat, tpw_data, 
                   cmap='viridis', shading='auto', alpha=0.8)
cb = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02, aspect=30)
cb.set_label('Total Precipitable Water (mm)', color='white', size=10)
cb.ax.yaxis.set_tick_params(color='white', labelcolor='white')

# Layer 2: RAP Advection Contours
# Positive values = Moisture Advection (Red)
clevs = np.arange(2, 22, 4)
cs = ax.contour(lons_r, lats_r, adv, levels=clevs, colors='red', linewidths=1.5)
ax.clabel(cs, inline=True, fontsize=10, fmt='%d', colors='white')

# Layer 3: RAP Wind Vectors
ax.quiver(lons_r[::3], lats_r[::3], u[::3, ::3], v[::3, ::3], 
          color='white', scale=450, width=0.003)

# Geography
ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='white', linewidth=0.8)
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), edgecolor='cyan')

plt.title(f"MIMIC-TPW2 + RAP 925mb Moisture Advection\nValid: {datetime.datetime.utcnow().strftime('%H:%MZ %d %b %Y')}", 
          color='white', fontsize=14, pad=15)

plt.savefig('output_map.png', facecolor='black', edgecolor='none', bbox_inches='tight', dpi=150)

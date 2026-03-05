import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
from siphon.catalog import TDSCatalog
import datetime
import numpy as np

# --- 1. SETUP TIME ---
# MIMIC files are usually valid on the top of the hour. 
# We look for the most recent available file.
now = datetime.datetime.utcnow()
current_hour = now.replace(minute=0, second=0, microsecond=0)

# --- 2. FETCH MIMIC-TPW2 (NC Sector) ---
# Trying the most common hourly path structure
try:
    # Pattern: mtpw2/data/YYYY/YYYYMMDD.nc (example structure)
    # If a direct NC file isn't available, we pull the global and crop it.
    tpw_url = f"https://tropic.ssec.wisc.edu/real-time/mtpw2/webAnims/tpw_nrl_colors/conus/current.nc"
    ds_tpw = xr.open_dataset(tpw_url)
    tpw_data = ds_tpw['tpw'].sel(lat=slice(28, 42), lon=slice(-90, -70))
except Exception as e:
    print(f"Primary Satellite Path Failed: {e}")
    # FALLBACK: Use the direct CONUS GIF for the background if netCDF is blocked
    # (This is much more stable than the NC link)
    tpw_data = None

# --- 3. FETCH RAP 13KM (925MB) ---
try:
    # Use 'isobaric' to match RAP coordinate naming
    cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml')
    ds_rap = cat.datasets[0].remote_access(use_xarray=True)
    subset_rap = ds_rap.metpy.sel(isobaric=925 * units.hPa, 
                                 lat=slice(42, 28), lon=slice(360-90, 360-70))
    
    u = subset_rap['u-component_of_wind_isobaric'].metpy.unit_array
    v = subset_rap['v-component_of_wind_isobaric'].metpy.unit_array
    q = subset_rap['Specific_humidity_isobaric'].metpy.unit_array
    lats_r, lons_r = subset_rap.lat, subset_rap.lon
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lons_r, lats_r)
    adv = mpcalc.advection(q, u, v, dx=dx, dy=dy) * 1e6 
except Exception as e:
    print(f"RAP Access Error: {e}")
    adv = None

# --- 4. PLOT ---
fig = plt.figure(figsize=(12, 9), facecolor='black')
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_extent([-88, -74, 31, 40])

if tpw_data is not None:
    im = ax.pcolormesh(tpw_data.lon, tpw_data.lat, tpw_data, cmap='viridis', alpha=0.8)
    plt.colorbar(im, ax=ax, label='TPW (mm)', pad=0.02)

if adv is not None:
    cs = ax.contour(lons_r, lats_r, adv, levels=np.arange(2, 22, 4), colors='red')
    ax.clabel(cs, inline=True, fontsize=10, fmt='%d')
    ax.quiver(lons_r[::3], lats_r[::3], u[::3, ::3], v[::3, ::3], color='white', scale=400)

ax.add_feature(cfeature.STATES, edgecolor='white')
plt.title(f"NC Regional Analysis\n{now.strftime('%H:%MZ %d %b %Y')}", color='white')
plt.savefig('output_map.png', facecolor='black', bbox_inches='tight')

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from copy import deepcopy
import rioxarray as rxr
import xarray as xr
import seaborn as sns
import colorsys
import cmocean
import xycmap
import dask
import os, sys, re

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf


# and whether or not to use the model that included the geo-coord polynomial
# main variable
var = 'NIRv'
# asynch neighborhood radius
neigh_rad = 100
# whether or not to use results of the model that included geocoods as covars
include_coords = 'y'
# whether or not to only use the top covars in the summary map
only_top_covars = False
# asynch percentile threshold below which to drop pixels from summary map
asynch_thresh = 75
# process the data raw, or read in processed files?
process_raw = False

# set path
data_dir = phf.EXTERNAL_RF_DATA_DIR

# load countries and subnational boundaries
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR, 'NewWorldFile_2020.shp'))
# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR,f)).to_crs(8857))
subnational = pd.concat(subnational)

# top covariates (in case only mapping those)
top_covars = ['ppt.asy', 'tmp.min.asy', 'veg.ent']

covar_cbar_labels = {
    'ppt.asy': '$\Delta dist_{seas_{P}}/\Delta  dist_{geo}$',
    'tmp.min.asy': '$\Delta dist_{seas_{T_{min}}}/\Delta  dist_{geo}$',
    'veg.ent': '$entropy$',
    'tmp.max.asy': '$\Delta dist_{seas_{T_{max}}}/\Delta  dist_{geo}$',
    'cld.asy': '$\Delta dist_{seas_{cld}}/\Delta  dist_{geo}$',
    'def.asy': '$\Delta dist_{seas_{def}}/\Delta  dist_{geo}$',
    'vrm.med': '$med_{VRM}$',
}

# create SHAP-interpretation map
# min-max scaler
def minmax_scale(vals, bot_pctile, top_pctile):
    scaled = (vals - np.nanpercentile(vals, bot_pctile))/(
                     np.nanpercentile(vals, top_pctile) -
                     np.nanpercentile(vals, bot_pctile))
    # clip to [0,1]
    scaled = np.clip(scaled, 0, 1)
    return scaled


# stack all covars' SHAP maps
files = [f for f in os.listdir(data_dir) if re.search(
    f'SHAP_map_{include_coords}COORDS_[cdptv].*_{var}_{neigh_rad}km.tif', f)]
# reorder with top covars first, so that map colors of top covars are same
# color regardless of whether other colors are included in the map
reordered_files = []
for covar in covar_cbar_labels.keys():
    for f in files:
        if re.search(f'{covar}_{var}', f):
            reordered_files.append(f)
files = reordered_files
if only_top_covars:
    files = [f for f in files if re.search(f'(?<=COORDS_).*(?=_{var})', f).group() in top_covars]
shap_maps = []
covars = []
for f in files:
    covar = re.search(f'(?<=COORDS_).*(?=_{var})', f).group()
    covars.append(covar)
    shap_map = rxr.open_rasterio(os.path.join(data_dir, f), masked=True)[0]
    shap_map.name = covar
    shap_maps.append(shap_map)
xr_stack = xr.concat(shap_maps, dim=pd.Index(covars, name='covar'))
da = dask.array.stack(xr_stack)

# rechunk (to run slowly but successfully on laptop)
# NOTE: axis sizes apparently needn't be exact multiples of their chunk sizes!
da = da.rechunk(chunks=(da.shape[0], 1, 53, 500))

if process_raw:
    # saturation: scaled standard deviation (as a metric of top covar predominance)
    # NOTE: 1 = max predominance; 0 = max evenness of importance
    std = deepcopy(shap_maps[0])
    std.values = minmax_scale(np.std(da, axis=0), 1, 99)
    assert np.nanmin(std) == 0 and np.nanmax(std) == 1
    # hue: predominant var
    predom = deepcopy(shap_maps[0])
    predom.values = np.argmax(np.abs(da), axis=0).compute()
    # NOTE: expressed in 360-degree values for hue (but cast as decimals,
    #       as required by the colorsys.hsv_to_rgb function I'm using), 
    #       yellow = 59/359, cyan = 179/359, and magenta = 299/359
    #predom = (59 + (120*predom))/359
    # mask to the std values
    predom = predom.where(pd.notnull(std))
    assert np.all(np.unique(predom)[~np.isnan(np.unique(predom))] ==
                  np.array([i/359 for i in [59,179,299]]))
else:
    if only_top_covars:
        file_suffix = 'top'
    else:
        file_suffix = 'all'
    std = rxr.open_rasterio(os.path.join(data_dir,
                    f'SHAP_std_{file_suffix}.tif'), masked=True)[0]
    predom = rxr.open_rasterio(os.path.join(data_dir,
                    f'SHAP_predom_{file_suffix}.tif'), masked=True)[0]

# predominance colormap
# NOTE: adapted from 'light' colormap at https://personal.sron.nl/~pault/#sec:qualitative
colors = [
          '#77AADD', # light blue
          '#FFAABB', # pink
          '#BBCC33', # pear
          '#EE8866', # orange
          '#EEDD88', # light yellow
          '#99DDFF', # light cyan
          '#44BB99', # mint
          #'#AAAA00', # olive
          #'#ADADAD', # light grey
         ]
cmap = ListedColormap(colors)

# value: asynchrony (will make lower-asynchrony pixels appear less pronounced
# because they'll be closer to black)
asynch = minmax_scale(rxr.open_rasterio(os.path.join(data_dir,
        f'{var}_STRICT_asynch_{neigh_rad}km.tif'), masked=True)[0], 0, 99)
# mask out pixels less than the percentile threshold
asynch = asynch.where(asynch >= np.nanpercentile(asynch, asynch_thresh))
# reproject and rescale
asynch = asynch.rio.reproject_match(std)
asynch = minmax_scale(asynch, 1, 99)
assert np.nanmin(asynch) == 0 and np.nanmax(asynch) == 1

# mask everything to the thresholded asynch map
std = std.where(pd.notnull(asynch))
predom = predom.where(pd.notnull(asynch))

# save rasters
if process_raw:
    if only_top_covars:
        file_suffix = 'top'
    else:
        file_suffix = 'all'
    for lyr_name, lyr in {'std': std, 'predom': predom}.items():
        lyr.rio.to_raster(os.path.join(data_dir,
                                       f'SHAP_{lyr_name}_{file_suffix}.tif'))


fig = plt.figure()
ax = fig.add_subplot()
divider = make_axes_locatable(ax)
cbar_ax = divider.append_axes('bottom', size='7%', pad=0.2)
predom.plot.imshow(ax=ax,
                   cmap=cmap,
                   add_colorbar=True,
                   cbar_ax=cbar_ax,
                   cbar_kwargs={'orientation': 'horizontal'},
                   zorder=1,
                  )
cbar_sec_width = (len(covars)-1)/len(covars)
ax.set_xticks(())
ax.set_yticks(())
cbar_ticks = np.linspace(cbar_sec_width/2,
                         (len(covars)-1)-(cbar_sec_width/2),
                         len(covars))
cbar_ax.set_xticks(cbar_ticks, covars)
for i in range(len(covars)):
    tmp_cmap = LinearSegmentedColormap.from_list('', ['#ffffff', colors[i]])
    std.where(predom==i).plot.imshow(ax=ax,
                                     cmap=tmp_cmap,
                                     add_colorbar=False,
                                     #vmin=0,
                                     #vmax=1,
                                     zorder=1+i
                                    )
subnational.to_crs(asynch.rio.crs).plot(ax=ax,
                                        color='none',
                                        edgecolor='black',
                                        linewidth=0.4,
                                        alpha=0.6,
                                        zorder=100,
                                       )
countries.to_crs(asynch.rio.crs).plot(ax=ax,
                                      color='black',
                                      edgecolor='black',
                                      linewidth=0,
                                      alpha=0.8,
                                      zorder=0,
                                     )
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_xlim(asynch.rio.bounds()[0::2])
ax.set_ylim(asynch.rio.bounds()[1::2])

fig.subplots_adjust(left=0.02,
                    right = 0.98,
                    bottom=0.08,
                    top=0.98,
                   )
fig.savefig('SHAP_summary_map.png', dpi=600)

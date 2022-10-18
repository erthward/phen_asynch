import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from copy import deepcopy
import rioxarray as rxr
import seaborn as sns
import colorsys
import cmocean
import dask
import os, sys, re

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf


# indicate the main variable and neighborhood radius,
# and whether or not to use the model that included the geo-coord polynomial
var = 'NIRv'
neigh_rad = 100
include_coords = 'y'

# set path
data_dir = phf.EXTERNAL_RF_DATA_DIR

top_covars = {'ppt.asy': f'pr_asynch_{neigh_rad}km.tif',
              'tmp.min.asy': f'tmmn_asynch_{neigh_rad}km.tif',
              'veg.ent': 'MODIS_IGBP_veg_entropy.tif',
             }
top_covar_cbar_labels = {
    'ppt.asy': '$\Delta dist_{seas_{P}}/\Delta  dist_{geo}$',
    'tmp.min.asy': '$\Delta dist_{seas_{T_{min}}}/\Delta  dist_{geo}$',
    'veg.ent': '$entropy$',
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


# stack all 3 top covars' SHAP maps
shap_maps = [rxr.open_rasterio(os.path.join(phf.data_dir,
        f'SHAP_map_{include_coords}COORDS_{tc}_{var}_{neigh_rad}km.tif'),
                               masked=True)[0] for tc in top_covars]
da = dask.array.stack(shap_maps)

# saturation: 1 - scaled standard deviation (as a metric top covar predominance)
std = deepcopy(shap_maps[0])
std.values = 1 - minmax_scale(np.std(da, axis=0), 1, 99)
assert np.nanmin(std) == 0 and np.nanmax(std) == 1

# hue: predominant var
predom = deepcopy(shap_maps[0])
predom.values = np.argmax(np.abs(da), axis=0).compute()
# NOTE: expressed in 360-degree values for hue (but cast as decimals,
#       as required by the colorsys.hsv_to_rgb function I'm using), 
#       yellow = 59/359, cyan = 179/359, and magenta = 299/359
predom = (59 + (120*predom))/359
# mask to the std values
predom = predom.where(pd.notnull(std))
assert np.all(np.unique(predom)[~numpy.isnan(np.unique(predom))] ==
                                    np.array([i/359 for i in [59,179,299]]))

# value: asynchrony (will make lower-asynchrony pixels appear less pronounced
# because they'll be closer to black)
asynch = minmax_scale(rxr.open_rasterio(os.path.join(phf.data_dir,
                        'NIRv_STRICT_asynch_100km.tif'), masked=True)[0], 0, 99)
asynch = asynch.rio.reproject_match(std)
asynch = minmax_scale(asynch, 1, 99)
assert np.nanmin(std) == 0 and np.nanmax(std) == 1

# mask everything to the asynch map
std = std.where(pd.notnull(asynch))
predom = predom.where(pd.notnull(asynch))

# stack all three, then convert to RGB
def hsv_ax_vals_2_rgb(vals):
    if np.sum(np.isnan(vals)) > 0:
        return np.array([np.nan]*3)
    else:
        return colorsys.hsv_to_rgb(h=vals[0], s=vals[1], v=vals[2])


hsv = dask.array.stack([predom, std, asynch]).rechunk([3,200,200])

rgb = xr.concat([deepcopy(std.expand_dims('color',
                                          0)) for _ in range(3)], dim='color')

# NOTE: GETTING A CONFUSING ERROR TRYING TO USE dask.array.apply_along_axis...
rgb.values = np.apply_along_axis(hsv_ax_vals_2_rgb, 0, hsv.compute())

# save to GeoTIFF
var_order = '_'.join([*top_covars])
rgb.rio.to_raster(os.path.join(phf.EXTERNAL_RF_DATA_DIR,
            f'SHAP_hsv_map_{var_order}_{include_coords}COORDS_{var}_{neigh_rad}km.tif'))


import rioxarray as rxr
import dask
import numpy as np
import pandas as pd
import geopandas as gpd
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os, sys, re

'''
TODO:
    - use asynch-clustering from final analysis to mask map instead of
      pixel-wise asynch values?
'''


# include lat and lon (i.e., y and x) maps?
include_latlon = False

# only use top-importance covars?
only_top_import_covars = False

# set top asynch percentile below which to mask out pixels
asynch_pctile = 90

# get var (SIF or NIRv) and neigh_rad (in km) for which to compare shap maps
var = sys.argv[1]
neigh_rad = int(sys.argv[2])

# get all valid files
data_dir = '/media/deth/SLAB/diss/3-phn/corr_data/'
shap_files = [f for f in os.listdir(data_dir) if f.endswith('.tif')]
shap_files = [f for f in shap_files if (var in f) and (str(neigh_rad)+'km' in f)]
if not include_latlon:
    shap_files = [f for f in shap_files if
                  not re.search('(?<=^SHAP_map_)[xy].?1?_%s' % var, f)]
if only_top_import_covars:
    shap_files = [f for f in shap_files if ('ppt' in f or
                                            'tmp.min' in f or
                            # NOTE: need x and y in case include_latlon == True
                            re.search('(?<=^SHAP_map_)[xy].?1?_%s' % var, f))]

# read all covars
covar_names = [re.search('(?<=^SHAP_map_).*(?=_%s)' % var,
                         f).group() for f in shap_files]
das = [rxr.open_rasterio(os.path.join(data_dir, f),
                         masked=True)[0] for f in shap_files]
covars = dask.array.stack(das)

# rechunk to make computable on my laptop
# NOTE: need factors of array sized (9, 2223, 6667)
covars = covars.rechunk(chunks=(9,741,113))


# collapse into single map with integer indicating layer with largest absolute
# SHAP value
# NOTE: any way to indicate SHAP val was + or -?
max_abs_shap = np.argmax(np.abs(covars), axis=0).compute()

# put back into a rioxarray object
out_da = deepcopy(das[0])
out_da[:,:] = max_abs_shap

# read in asynch file and mask to match its footprint
asynch_data_dir = '/media/deth/SLAB/diss/3-phn/GEE_outputs/final/'
asynch = rxr.open_rasterio(os.path.join(asynch_data_dir,
                    '%s_asynch_%ikm.tif' % (var, neigh_rad)), masked=True)[0]

# TODO: DELETE ME
asynch = asynch.rio.reproject_match(out_da)

out_da = out_da.where(pd.notnull(asynch), np.nan)

# mask to only top Nth percentile of asynch values
out_da = out_da.where(asynch>=np.nanpercentile(asynch, asynch_pctile))

# set colormap for the covars
# NOTE: using 12-color palette at http://tsitsul.in/blog/coloropt/
color_hex = ['#ebac23',
             '#b80058',
             '#008cf9',
             '#006e00',
             '#00bbad',
             '#d163e6',
             '#b24502',
             '#ff9287',
             '#5954d6',
             '#00c6f8',
             '#878500',
             '#00a76c',
             '#bdbdbd',
            ]
cmap = ListedColormap(color_hex[:len(covar_names)])
# make missing values black (to make colors more discernable)
#cmap.set_bad('#ffffff')

# load countries boundaries
countries_dir = ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                 'seasonal_asynchrony/data/bounds')
countries = gpd.read_file(os.path.join(countries_dir, 'NewWorldFile_2020.shp'))

# plot it
fig = plt.figure(figsize=(22, 10))
ax = fig.add_subplot(111)
countries.to_crs(out_da.rio.crs).plot(color='black',
                                      ax=ax,
                                      zorder=0)
im = out_da.plot.imshow(cmap=cmap, ax=ax, zorder=1)
cbar = im.colorbar
ymin, ymax = cbar.ax.get_ylim()
tick_locs = np.linspace(ymax/len(covar_names)/2,
                        ymax - (ymax/len(covar_names)/2),
                        len(covar_names))
cbar.ax.set_yticks(tick_locs, covar_names)
# format
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('predomnant covars, by |SHAP| (%s, %i km)' % (var, neigh_rad))
fig.subplots_adjust(left=0.01, right=1.06, bottom=0.01, top=0.95)
fig.savefig('SHAP_interp_map_%s_%ikm%s%s.png' % (var, neigh_rad,
                            '_LATLON' * include_latlon,
                            '_TOPIMPORT' * only_top_import_covars), dpi=700)

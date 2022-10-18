import rioxarray as rxr
import dask
import numpy as np
import pandas as pd
import geopandas as gpd
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os, sys, re


# include lat and lon (i.e., y and x) maps?
include_latlon = False

# only use top-importance covars?
use_only_top_import_covars = True

# set top asynch percentile below which to mask out pixels
asynch_pctile = 75

# get var (SIF or NIRv) and neigh_rad (in km) for which to compare shap maps
var = sys.argv[1]
neigh_rad = int(sys.argv[2])
include_coords = sys.argv[3]

# get all valid files
data_dir = '/media/deth/SLAB/diss/3-phn/rf_data/'
shap_files = [f for f in os.listdir(data_dir) if (f.startswith('SHAP_map')
                                                  and f.endswith('.tif'))]
shap_files = [f for f in shap_files if (var in f and
                                        str(neigh_rad)+'km' in f and
                                        f'{include_coords}COORDS' in f)]
if not include_latlon:
    shap_files = [f for f in shap_files if not
                      re.search('(?<=^SHAP_map_[yn]COORDS_)[xy].?1?_%s' % var, f)]
if use_only_top_import_covars:
    shap_files = [f for f in shap_files if ('ppt' in f or
                                            'tmp.min' in f or
                            # NOTE: need x and y in case include_latlon == True
                            re.search('(?<=^SHAP_map_[yn]COORDS_)[xy].?1?_%s' % var, f))]

# sort the shap_files list (to ensure identical order every time I run this
# and then the follow-on figure-producing script
shap_files = np.sort(shap_files)

# read all covars
covar_names = [re.search('(?<=^SHAP_map_[yn]COORDS_).*(?=_%s)' % var,
                         f).group() for f in shap_files]
das = [rxr.open_rasterio(os.path.join(data_dir, f),
                         masked=True)[0] for f in shap_files]
covars = dask.array.stack(das)

# rechunk to make computable on my laptop
# NOTE: need factors of array sized (9, 2223, 6667)
covars = covars.rechunk(chunks=(9,741,113))

# make sure each cell is either all valid values or all NAs
na_sums = np.sum(np.isnan(covars), axis=0).compute()
assert [*np.unique(na_sums)] == [0, len(covar_names)]

# collapse into single map with integer indicating layer with largest absolute
# SHAP value
# NOTE: any way to indicate SHAP val was + or -?
max_abs_shap = np.argmax(np.abs(covars), axis=0).compute()

# put back into a rioxarray object
out_da = deepcopy(das[0])
out_da[:,:] = max_abs_shap

# mask out places with NaNs
out_da.values[na_sums>0] = np.nan

assert np.nanmin(out_da) == 0
assert np.nanmax(out_da) == (len(covar_names)-1)

# read in asynch file and mask to match its footprint
asynch_data_dir = '/media/deth/SLAB/diss/3-phn/GEE_outputs/final/'
asynch = rxr.open_rasterio(os.path.join(asynch_data_dir,
                    '%s_asynch_%ikm.tif' % (var, neigh_rad)), masked=True)[0]

# TODO: DELETE ME
asynch = asynch.rio.reproject_match(out_da)

out_da = out_da.where(pd.notnull(asynch), np.nan)

# mask to only top Nth percentile of asynch values
out_da = out_da.where(asynch>=np.nanpercentile(asynch, asynch_pctile))

# write to file
filename = 'SHAP_interp_map_%s_%ikm%s%s.tif' % (var, neigh_rad,
                                                '_LATLON' * include_latlon,
                                        '_TOPIMPORT' * use_only_top_import_covars)
out_da.rio.to_raster(os.path.join(data_dir, filename))




# plot map for cross-checking with final fig
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

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
im = out_da.plot.imshow(cmap=cmap,
                        ax=ax,
                        add_colorbar=True,
                        vmin=0,
                        vmax=len(covar_names)-1,
                       )
cbar = im.colorbar
axmin, axmax = cbar.ax.get_ylim()
tick_locs = np.linspace(axmax/len(covar_names)/2,
                        axmax - (axmax/len(covar_names)/2),
                        len(covar_names))
cbar.ax.set_yticks(tick_locs, covar_names)
# format
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('')
cbar.ax.tick_params(length=0, labelsize=12)

for tick in cbar.ax.get_yticklabels():
    tick.set_rotation(0)

# write cross-check figure to file
filename = 'SHAP_interp_map_%s_%ikm%s%s_CROSSCHECK.png' % (var, neigh_rad,
                                                '_LATLON' * include_latlon,
                                        '_TOPIMPORT' * use_only_top_import_covars)
fig.savefig(os.path.join(data_dir, filename), dpi=700)



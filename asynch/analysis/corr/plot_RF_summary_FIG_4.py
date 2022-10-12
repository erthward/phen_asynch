import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from copy import deepcopy
import rioxarray as rxr
import seaborn as sns
import dask
import os, sys, re

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf



"""
TODO:
    - get running successfully after re-rasterizing SHAP vals
    - clip all maps to the same bounding box used in other figs
    - save, then add, other importance metric
    - add letters to label parts
    - other TODOs
    - if keeping interp_map, could like colors in the map to font colors or
    something used to label rows at the right
    - in interp map, use asynch-clustering from final analysis to mask map instead of
      pixel-wise asynch values?
"""


# set plotting params
title_fontdict = {'fontsize': 36}
axlabel_fontdict = {'fontsize': 44}
ticklabel_fontsize=20

# set paths
data_dir = phf.EXTERNAL_CORR_DATA_DIR

# load shapefiles
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                       'NewWorldFile_2020.shp'))
# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR,f)).to_crs(8857))
subnational = pd.concat(subnational)

var = 'NIRv'
neigh_rad = 100
# TODO: var filenames
top_covars = {'tmp.min.nsd': ('CHELSA_bio6_1981-2010_V.2.1_5km_'
                              '10CELLRAD_NEIGHSD.tif'),
              'ppt.sea.nsd': ('CHELSA_bio15_1981-2010_V.2.1_5km_'
                              '10CELLRAD_NEIGHSD.tif'),
              'ppt.asy': f'pr_asynch_{neigh_rad}km.tif',
              'tmp.min.asy': f'tmmn_asynch_{neigh_rad}km.tif',
             }


# set up figure
fig = plt.figure(figsize=(28, 21))
gs = fig.add_gridspec(ncols=120, nrows=90)


#------------------------------------------
# 1. plot original vars and their SHAP maps
row_idxs = np.arange(5, 90, 20)
col_idxs = np.arange(80, 121, 20)
for i, top_covar_item in enumerate(top_covars.items()):
    top_covar, covar_filename = top_covar_item
    shap_filename = f'SHAP_map_{top_covar}_{var}_{neigh_rad}km.tif'
    ax_covar = fig.add_subplot(gs[row_idxs[i]:row_idxs[i+1],
                                  col_idxs[0]:col_idxs[1]])
    ax_shap = fig.add_subplot(gs[row_idxs[i]:row_idxs[i+1],
                                 col_idxs[1]:col_idxs[2]])
    j = 0
    for ax, fn in zip([ax_covar, ax_shap], [covar_filename, shap_filename]):
        divider = make_axes_locatable(ax)
        where = 'right'
        orientation = 'vertical'
        size = '4%'
        cax = divider.append_axes(where, size=size, pad=0.2)
        rast = rxr.open_rasterio(os.path.join(data_dir, fn), masked=True)[0]
        rast = rast.rio.write_crs(4326).rio.reproject(8857)
        # NOTE: hack around the stupid xarray AttributeError
        rast.attrs['long_name'] = ''
        rast.plot.imshow(ax=ax,
                         zorder=0,
                         cmap='cividis',
                         vmin=np.nanpercentile(rast, 1),
                         vmax=np.nanpercentile(rast, 99),
                         add_colorbar=True,
                         cbar_ax=cax,
                         cbar_kwargs = {'orientation': orientation},
                        )
        cax.tick_params(labelsize=20)
        cax.set_ylabel(var, fontdict={'fontsize': 26})
        subnational.to_crs(8857).plot(ax=ax,
                                      color='none',
                                      edgecolor='black',
                                      zorder=1,
                                      alpha=0.6,
                                     )
        countries.to_crs(8857).plot(ax=ax,
                                    color='none',
                                    edgecolor='black',
                                    linewidth=1,
                                    alpha=0.8,
                                    zorder=2,
                                   )
        # format axes
        ax.set_xlim(rast.rio.bounds()[0::2])
        ax.set_ylim(rast.rio.bounds()[1::2])
        # NOTE: chopping off western edge because the equal earth projection
        #       makes NZ appear twice
        ax.set_xlim(0.95 * ax.get_xlim()[0], ax.get_xlim()[1])
        if i == 0:
            if j == 0:
                ax.set_title('covariate')
            else:
                ax.set_title('SHAP values', fontdict=axlabel_fontdict)
        else:
            ax.set_title('')
        if j == 0:
            ax.set_ylabel(var, fontdict=axlabel_fontdict)
        else:
            ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_xticks(())
        ax.set_yticks(())
        del rast
        j += 1


#----------------------------------
# 2. plot covar importance bar plot

dfs = []
# TODO: ADD IN OTHER IMPORTANCE METRIC
#for import_metric in ['SHAP', 'info']:
for import_metric in ['SHAP', 'SHAP']:
    import_df = pd.read_csv(os.path.join(data_dir,
                    f'rf_{import_metric}_importance_{var}_{neigh_rad}km.csv'))
    import_df.columns = ['covariate', 'importance']
    import_df.sort_values(by='importance', ascending=False, inplace=True)
    dfs.append(import_df)

col_idx = 0
for i, df in enumerate(dfs):
    ax = fig.add_subplot(gs[:30, col_idx:(col_idx+30)])
    sns.barplot(y='covariate',
                x='importance',
                data=df,
                hue='importance',
                orient='h',
                palette='OrRd',
                ax=ax,
               )
    ax.set_xlabel('RF covariate', fontdict=axlabel_fontdict)
    ax.set_ylabel('importance', fontdict=axlabel_fontdict)
    ax.set_title('')
    ax.tick_params(labelsize=ticklabel_fontsize)
    # TODO:
    # add labels for parts A. and B.
    #fig.axes[0].text(-5.8, -5, 'A.', size=24, weight='bold')
    #fig.axes[-2].text(1.11*fig.axes[-2].get_xlim()[0],
    #                 1.065*fig.axes[-2].get_ylim()[1],
    #                 'B.', size=24, weight='bold')
    col_idx += 30


#--------------------------------
# 3. plot SHAP interpretation map

ax = fig.add_subplot(gs[30:, :80])
divider = make_axes_locatable(ax)
where = 'bottom'
orientation = 'vertical'
size = '7%'
cax = divider.append_axes(where, size=size, pad=0.2)

# include lat and lon (i.e., y and x) maps?
include_latlon = False

# only use top-importance covars?
use_only_top_import_covars = False

# use a special color for the background?
bg_color = None
#bg_color = '#ffffff'

# set top asynch percentile below which to mask out pixels
asynch_pctile = 90

# get all valid files
shap_files = [f for f in os.listdir(data_dir) if (f.startswith('SHAP_map') and
                                                  f.endswith('.tif'))]
shap_files = [f for f in shap_files if (var in f) and (str(neigh_rad)+'km' in f)]
if not include_latlon:
    shap_files = [f for f in shap_files if
                  not re.search('(?<=^SHAP_map_)[xy].?1?_%s' % var, f)]
if use_only_top_import_covars:
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
asynch_data_dir = phf.EXTERNAL_DATA_DIR
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
if bg_color is not None:
    cmap.set_bad(bg_color)

# plot it
countries.to_crs(out_da.rio.crs).plot(color='black',
                                      ax=ax,
                                      zorder=0)
im = out_da.plot.imshow(cmap=cmap,
                        ax=ax,
                        add_colorbar=True,
                        cbar_ax=cax,
                        zorder=1,
                       )
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
ax.set_title('predomnant covars, by |SHAP value|', fontdict=title_fontdict)

# adjust subplots and save
fig.subplots_adjust(bottom=0.03,
                    top=0.92,
                    left=0.08,
                    right=0.9,
                    wspace=0,
                    hspace=0,
                   )
fig.savefig('FIG_4_RF_summary.png', dpi=700)

"""
TODO:
    - clean up cruft
    - uncomment and figure out lower code for viz and RF supp figs
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from matplotlib.patches import Polygon as mplPolygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as path_effects
from collections import Counter as C
from palettable.cartocolors.sequential import PinkYl_7
from palettable.scientific.sequential import Acton_20_r
from shapely.geometry import Polygon as shapelyPolygon
from math import pi
from copy import deepcopy
import seaborn as sns
import statsmodels.api as sm
import palettable
import cmocean
import colorsys
import cmocean
import dask
import sys
import re
import os

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf


# set data-directory path
data_dir = phf.EXTERNAL_RF_DATA_DIR

# general plotting params:
title_fontsize = 12
rowlab_fontsize = 18
axislab_fontsize = 13
ticklab_fontsize = 12
cbarlab_fontsize = 12
cbar_axlab_fontsize = 35
cbar_ticklab_fontsize = 24
annot_fontsize = 14
scat_label_fontsize = 10
scat_label_fontcolor = 'black'
scat_label_linewidth = 3
fig_width = 10.5
fig_height = 12
dpi = 400
n_ticklabels = 5
subplots_adj_left=0.06
subplots_adj_bottom=0.05
subplots_adj_right=0.96
subplots_adj_top=0.95
subplots_adj_wspace=0.14
subplots_adj_hspace=0.30
central_curve_color = '#220000'
central_linewidth = 2
central_alpha = 1
neighbor_linewidth = 0.1
neighbor_alpha = 0.5
rad_linestyle = ':'
rad_color = 'white'
rad_linewidth = 2
rad_alpha = 0.75
rad_mask_alpha = 0.175
supp_axlabel_fontdict = {'fontsize': 18}
cbarlabel_fontdict = {'fontsize': 14}
cbarticklabel_fontsize = 12
partlabel_fontsize=40

# TODO: FIX EQUAL EARTH PROJECTION?
plot_crs = 8857


# asynch params for main results:
# main variable
var = 'NIRv'
# asynch neighborhood radius
neigh_rad = 100
# whether or not to use results of the model that included geocoods as covars
include_coords = 'y'

# SHAP-plotting params:
# whether or not to only use the top covars in the summary map
only_top_covars = True
# asynch percentile threshold below which to drop pixels from summary map
asynch_thresh = 75
# process the data raw, or read in processed files?
process_raw = False

# asynch colormap
asynch_cmap = 'cmo.dense_r'

# and predominance colormap
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
if only_top_covars:
    colors = colors[:3]
predom_cmap = ListedColormap(colors)


# set up SHAP data:
def minmax_scale(vals, bot_pctile, top_pctile):
    """
    min-max scaler
    """
    scaled = (vals - np.nanpercentile(vals, bot_pctile))/(
                     np.nanpercentile(vals, top_pctile) -
                     np.nanpercentile(vals, bot_pctile))
    # clip to [0,1]
    scaled = np.clip(scaled, 0, 1)
    return scaled


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

covar_longnames = {
    'ppt.asy': 'precipitation\nasynchrony',
    'tmp.min.asy': 'min. temperature\nasynchrony',
    'veg.ent': 'vegetation\nentropy',
    'tmp.max.asy': 'max. temperature\nasynchrony',
    'cld.asy': 'cloud\nasynchrony',
    'def.asy': 'CWD\nasynchrony',
    'vrm.med': 'median vector\nruggedness metric',
    'no.predom': 'no predominant\ncovariate',
}

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
all_files = reordered_files
if only_top_covars:
    files = [f for f in all_files if re.search(f'(?<=COORDS_).*(?=_{var})', f).group() in top_covars]
else:
    files = all_files

def map_asynch(fig, cbar_axlab,
               gs=None, main_fig=True, var='NIRv',
               cbar_axlab_fontsize=cbar_axlab_fontsize,
               cbar_ticklab_fontsize=cbar_ticklab_fontsize):

    assert var in ['NIRv', 'SIF', 'tmmn', 'tmmx', 'pr', 'def', 'cloud']

    if var in ['NIRv', 'SIF']:
        files = [f for f in os.listdir(phf.EXTERNAL_DATA_DIR) if
                                        re.search('%s_STRICT_asynch' % var, f)]
    else:
        files = [f for f in os.listdir(phf.EXTERNAL_DATA_DIR) if
                                        re.search('%s_asynch' % var, f)]

    # cut down to just one file, if this is for the main fig
    if main_fig:
        files = [f for f in files if re.search('asynch_100km', f)]
    # otherwise arrange in top-down order of increasing neighborhood radius
    else:
        reordered_files = []
        for neigh_rad in [50, 100, 150]:
            neigh_rad_file = [f for f in files if re.search('asynch_%ikm' %
                                                            neigh_rad, f)]
            assert len(neigh_rad_file) == 1
            reordered_files.append(neigh_rad_file[0])
        files = reordered_files

    for ax_ct, file in enumerate(files):
        # get the neighborhood radius
        neigh_rad = int(re.search('(?<=asynch_)\d{2,3}(?=km\.tif)',
                                  file).group())

        # either grab the upper half of the main fig
        if main_fig:
            ax = fig.add_subplot(gs[:130, :])
        # or grab the next row of the supp fig
        else:
            ax = fig.add_subplot(3,1,ax_ct+1)

        # partition off a separate axis for the colormap
        divider = make_axes_locatable(ax)
        if main_fig:
            where = 'bottom'
            orientation = 'horizontal'
            size = '7%'
        else:
            where = 'right'
            orientation = 'vertical'
            size = '4%'
        cax = divider.append_axes(where, size=size, pad=0.2)

        # add jurisdictional boundaries
        if main_fig:
            plot_juris_bounds(ax, 0, 2, 3)

        # read in the raster data and prepare it
        rast = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                              file), masked=True)[0]
        rast = rast.rio.write_crs(4326).rio.reproject(plot_crs, nodata=np.nan)
                # NOTE: annoying AttributeError is because da.attrs['long_name']
        #       is retained as a tuple of names (rather than being subsetted
        #       by indexing) when I index a single layer out of an
        #       xarray.core.dataarray.DataArray;
        #       for now, a hacky fix is just assigning a string to that attr
        rast.attrs['long_name'] = ''
        rast.plot.imshow(ax=ax,
                         cmap=asynch_cmap,
                         vmin=np.nanpercentile(rast, 1),
                         vmax=np.nanpercentile(rast, 99),
                         add_colorbar=True,
                         cbar_ax=cax,
                         cbar_kwargs = {'orientation': orientation},
                         zorder=1,
                        )
        cax.tick_params(labelsize=cbar_ticklab_fontsize)
        if main_fig:
            cax.set_xlabel(cbar_axlab, fontdict={'fontsize': cbar_axlab_fontsize})
        else:
            cax.set_ylabel(cbar_axlab, fontdict={'fontsize': cbar_axlab_fontsize})
                # format axes
        ax.set_xlim(rast.rio.bounds()[0::2])
        ax.set_ylim(rast.rio.bounds()[1::2])
        # NOTE: chopping off western edge because the equal earth projection
        #       makes NZ appear twice
        ax.set_xlim(0.95 * ax.get_xlim()[0], ax.get_xlim()[1])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks(())
        ax.set_yticks(())
        # add axis title, if not main figure
        if main_fig:
            ax.set_title('')
        else:
            ax.set_title('%i km neighborhood' % neigh_rad,
                         fontdict={'fontsize': 21})

        del rast

if __name__ == '__main__':


#---------------------------------------
# plot original vars and their SHAP maps

    # get all covars
    all_covars = []
    shap_maps = []
    for f in all_files:
        covar = re.search(f'(?<=COORDS_).*(?=_NIRv)', f).group()
        all_covars.append(covar)
        shap_map = rxr.open_rasterio(os.path.join(data_dir, f), masked=True)[0]
        shap_map.name = covar
        shap_maps.append(shap_map)
    xr_stack = xr.concat(shap_maps, dim=pd.Index(all_covars, name='covar'))
    da = dask.array.stack(xr_stack)

    # loop over all SHAP tifs, to determine the min and max SHAP vals for plotting
    vmins = []
    vmaxs = []
    for j, covar in enumerate(all_covars):
        shap_filename = f'SHAP_map_{include_coords}COORDS_{covar}_{var}_{neigh_rad}km.tif'
        rast = rxr.open_rasterio(os.path.join(data_dir,
                                              shap_filename), masked=True)[0]
        vmin = np.nanpercentile(rast, 1)
        vmax = np.nanpercentile(rast, 99)
        vmins.append(vmin)
        vmaxs.append(vmax)
    vminmax = np.max(np.abs(vmins+vmaxs))
    shap_vmin = -vminmax
    shap_vmax = vminmax


    def change_cbar_ticklabels_to_scinot(cax, orientation):
        if orientation == 'horizontal':
            axis = 'x'
        else:
            axis = 'y'
        ticks = getattr(cax, f'get_{axis}ticks')()
        new_ticklabs = ['%0.1e' % t for t in ticks]
        getattr(cax, f'set_{axis}ticks')(ticks, new_ticklabs)


    fig_all_var_maps = plt.figure(figsize=(8,10.5))
    gs = fig_all_var_maps.add_gridspec(ncols=2, nrows=len(all_covars))

    for i, covar in enumerate(all_covars):
        covar_filename = all_files[i]
        shap_filename = f'SHAP_map_{include_coords}COORDS_{covar}_{var}_{neigh_rad}km.tif'
        ax_covar = fig_all_var_maps.add_subplot(gs[i, 0])
        ax_shap = fig_all_var_maps.add_subplot(gs[i, 1])
        # NOTE: cmo.dense_r matches the asynchrony map in fig 3
        cmaps = ['cmo.dense_r', 'cmo.curl_r']
        i = 0
        for ax, filename in zip([ax_covar, ax_shap],
                                [covar_filename, shap_filename]):
            divider = make_axes_locatable(ax)
            where = 'bottom'
            orientation = 'horizontal'
            size = '4%'
            cax = divider.append_axes(where, size=size, pad=0.2)
            rast = rxr.open_rasterio(os.path.join(data_dir,
                                                  filename), masked=True)[0]
            rast = rast.rio.write_crs(4326).rio.reproject(plot_crs,
                                                          nodata=np.nan)
            # NOTE: get rid of egregiously large values about longitudinal bounds
            #       of Equal Earth projection
            #rast = rast.where(np.abs(rast<=1e37))
            # NOTE: hack to get around the stupid xarray AttributeError
            rast.attrs['long_name'] = ''
            if i == 1:
                vmin = shap_vmin
                vmax = shap_vmax
            else:
                vmin = np.nanpercentile(rast, 1)
                vmax = np.nanpercentile(rast, 99)
            rast.plot.imshow(ax=ax,
                             zorder=0,
                             cmap=cmaps[i],
                             vmin=vmin,
                             vmax=vmax,
                             add_colorbar=True,
                             cbar_ax=cax,
                             cbar_kwargs = {'orientation': orientation},
                            )
            cax.tick_params(labelsize=cbarticklabel_fontsize)
            if i == 0:
                cbarlabel = covar
            else:
                cbarlabel = 'SHAP val'
            if orientation == 'vertical':
                axis = 'y'
            else:
                axis = 'x'
            getattr(cax, f'set_{axis}label')(cbarlabel, cbarlabel_fontdict)
            ax.set_xlim(rast.rio.bounds()[::2])
            ax.set_ylim(rast.rio.bounds()[1::2])
            if j == 0:
                if i == 0:
                    ax.set_ylabel('covariate', fontdict=supp_axlabel_fontdict)
                else:
                    ax.set_ylabel('influence', fontdict=supp_axlabel_fontdict)
            else:
                ax.set_ylabel('')
            if i == 0:
                lab = ax.set_title(covar, fontdict=supp_axlabel_fontdict)
            else:
                ax.set_title('')
            ax.set_xlabel('')
            ax.set_xticks(())
            ax.set_yticks(())
            del rast

            # add part label
            if i == 0 and j == 0:

                ax.text(1.35 * ax.get_xlim()[0],
                        1.4 * ax.get_ylim()[1],
                        'B.', size=partlabel_fontsize, weight='bold')

            elif i == 1 and j == 0:

                ax.text(1.35 * ax.get_xlim()[0],
                        1.4 * ax.get_ylim()[1],
                        'C.', size=partlabel_fontsize, weight='bold')

            # set number of x-axis ticks on cbar axes
            if i == 0:
                nticks = 4
            else:
                nticks = 5
            cax.xaxis.set_major_locator(plt.MaxNLocator(nticks))

            i += 1


    # adjust subplots and save
    fig_all_var_maps.subplots_adjust(bottom=0.02,
                        top=0.99,
                        left=0.06,
                        right=0.92,
                        wspace=0.02,
                        hspace=0.03,
                       )
    fig_all_var_maps.savefig('FIG_S15_all_var_maps.png', dpi=500)

    del fig_all_var_maps



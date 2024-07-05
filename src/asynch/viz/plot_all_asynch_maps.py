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
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf


# set data-directory path
data_dir = phf.EXTERNAL_RF_DATA_DIR

# general plotting params:
cbar_axlab_fontsize = 35
cbar_ticklab_fontsize = 24
dpi = 600
plot_crs = 8857

# vars to plot asynch for
vars = ['NIRv', 'SIF', 'tmmn', 'tmmx', 'pr', 'def', 'cloud']

# asynch colormap
asynch_cmap = palettable.scientific.sequential.LaJolla_20_r.mpl_colormap

cbar_axlab_dict = {'NIRv main': 'LSP asynchrony',
                   'NIRv': '$NIR_{V}\ asynch\ (\Delta NIR_{V}/\Delta m)$',
                   'SIF': '$SIF\  asynch\ (\Delta (mW m^{-2} sr^{-1} nm^{-1})/\Delta m)$',
                   'tmmn': '$tmp_{min}\ asynch\ (\Delta ^{\circ} C/\Delta m)$',
                   'tmmx': '$tmp_{max}\ asynch\ (\Delta ^{\circ} C/\Delta m)$',
                   'pr': '$ppt\ asynch\ (\Delta mm/\Delta m)$',
                   'def': '$cwd\ asynch\ (\Delta mm/\Delta m)$',
                   'cloud': '$cloud\ asynch\ (\Delta \%\ cover/\Delta m)$',

                  }


def plot_juris_bounds(ax, bg_zorder=0, subnat_zorder=1, nat_zorder=2,
                      polys_color='#060606',
                     ):
    """
    plot national and subnational jurisdictional bounds
    """
    if bg_zorder is not None:
        phf.plot_juris_bounds(ax,
                              lev1=False,
                              lev0_color=polys_color,
                              lev0_linewidth=0,
                              lev0_alpha=0.2,
                              lev0_zorder=bg_zorder,
                              crs=plot_crs,
                              strip_axes=False,
                             )
    phf.plot_juris_bounds(ax,
                          lev0_linecolor='gray',
                          lev0_linewidth=0.5,
                          lev0_alpha=0.8,
                          lev0_zorder=nat_zorder,
                          lev1_linecolor='gray',
                          lev1_linewidth=0.3,
                          lev1_alpha=0.7,
                          lev1_zorder=subnat_zorder,
                          crs=plot_crs,
                          strip_axes=False,
                         )


def map_asynch(fig, cbar_axlab,
               gs=None, var='NIRv',
               cbar_axlab_fontsize=cbar_axlab_fontsize,
               cbar_ticklab_fontsize=cbar_ticklab_fontsize):
    """
    plot the a phenological or climatic asynchrony map
    """
    assert var in ['NIRv', 'SIF', 'tmmn', 'tmmx', 'pr', 'def', 'cloud']
    if var in ['NIRv', 'SIF']:
        files = [f for f in os.listdir(phf.EXTERNAL_DATA_DIR) if
                                        re.search('%s_STRICT_asynch' % var, f)]
    else:
        files = [f for f in os.listdir(phf.EXTERNAL_DATA_DIR) if
                                        re.search('%s_asynch' % var, f)]
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
        # grab the next row of the supp fig
        ax = fig.add_subplot(3,1,ax_ct+1)
        # partition off a separate axis for the colormap
        divider = make_axes_locatable(ax)
        where = 'right'
        orientation = 'vertical'
        size = '4%'
        cax = divider.append_axes(where, size=size, pad=0.2)
        plot_juris_bounds(ax, 0, 2, 3)
        # read in the raster data and prepare it
        rast = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                              file), masked=True)[0]
        rast_proj = rast.rio.write_crs(4326).rio.reproject(plot_crs)
        rast = phf.mask_xarr_to_other_xarr_bbox(rast_proj, rast)
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
        cax.set_ylabel(cbar_axlab, fontdict={'fontsize': cbar_axlab_fontsize})
        # format axes
        ax.set_xlim(rast.rio.bounds()[0::2])
        ax.set_ylim(rast.rio.bounds()[1::2])
        phf.strip_axes_labels_and_ticks(ax)
        ax.set_title('%i km neighborhood' % neigh_rad,
                     fontdict={'fontsize': 21})
        del rast



# make both vars' supp figs (each one stacking all 3 neighborhood radii)
for n, var in enumerate(vars):
    print('\n\nNOW PRODUCING SUPPLEMENTAL FIG FOR %s..\n\n' % var)
    fig_supp = plt.figure(figsize=(9,12))
    map_asynch(fig_supp, cbar_axlab_dict[var],
               gs=None, var=var,
               cbar_axlab_fontsize=13, cbar_ticklab_fontsize=10)
    fig_supp.subplots_adjust(bottom=0.02, top=0.95, left=0.02, right=0.88)
    fig_supp.savefig(os.path.join(phf.FIGS_DIR, 'FIG_SUPP_%s_asynch_maps.png' %
                                  (var)), dpi=600)
    del fig_supp




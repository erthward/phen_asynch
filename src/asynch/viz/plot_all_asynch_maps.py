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
import matplotlib.ticker as ticker
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


# general plotting params:
cbar_axlab_fontsize = 7
cbar_ticklab_fontsize = 5
dpi = 600
plot_crs = 8857

# vars to plot asynch for
vars = ['NIRv', 'SIF', 'tmmn', 'tmmx', 'pr', 'def', 'cloud']

# asynch colormap
asynch_cmap = palettable.scientific.sequential.LaJolla_20_r.mpl_colormap

cbar_axlab_dict = {'NIRv main': 'LSP asynchrony',
                   'NIRv': ('$NIR_{V}\ asynch$'
                            '\n'
                            '$(\Delta NIR_{V}/\Delta m)$'),
                   'SIF': ('$SIF\ asynch$'
                           '\n'
                           '$(\Delta mW m^{-2} sr^{-1} nm^{-1}/\Delta m)$'),
                   'tmmn': ('$tmp_{min}\ asynch$'
                            '\n'
                            '$(\Delta ^{\circ} C/\Delta m)$'),
                   'tmmx': ('$tmp_{max}\ asynch$'
                            '\n'
                            '$(\Delta ^{\circ} C/\Delta m)$'),
                   'pr': ('$ppt\ asynch$'
                          '\n'
                          '$(\Delta mm/\Delta m)$'),
                   'def': ('$cwd\ asynch$'
                           '\n'
                           '$(\Delta mm/\Delta m)$'),
                   'cloud': ('$cloud\ asynch$'
                             '\n'
                             '$(\Delta \%\ cover/\Delta m)$'),
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
                          lev0_linewidth=0.1,
                          lev0_alpha=0.8,
                          lev0_zorder=nat_zorder,
                          lev1_linecolor='gray',
                          lev1_linewidth=0.1,
                          lev1_alpha=0.7,
                          lev1_zorder=subnat_zorder,
                          crs=plot_crs,
                          strip_axes=False,
                         )


def map_asynch(fig,
               gs,
               row,
               cbar_axlab,
               var='NIRv',
               cbar_axlab_fontsize=cbar_axlab_fontsize,
               cbar_ticklab_fontsize=cbar_ticklab_fontsize,
              ):
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
        ax = fig.add_subplot(gs[row, ax_ct])
        # formatter function for scientific notation on colorbars that need it
        def format_cbar(x, pos):
            a, b = '{:.1e}'.format(x).split('e')
            a = float(a)
            b = int(b)
            if a == 0:
                return "0"
            else:
                return "$%0.1f\\times10^{%i}$" % (a, b)
        # get cbar ax, if this is the last column
        if ax_ct == 2:
            cax = fig.add_subplot(gs[row, ax_ct+2])
            add_colorbar = True
            cbar_kwargs = {'orientation': 'vertical',
                          }
            # add formatter function for rows that need it
            if row in [1, 2, 3, 5]:
                cbar_kwargs['format'] = ticker.FuncFormatter(format_cbar)
        else:
            cax = None
            add_colorbar = False
            cbar_kwargs = None
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
                         #vmin=np.nanpercentile(rast, 1),
                         vmin=0,
                         vmax=np.nanpercentile(rast, 99),
                         add_colorbar=add_colorbar,
                         cbar_ax=cax,
                         cbar_kwargs=cbar_kwargs,
                         zorder=1,
                        )
        if cax is not None:
            cax.tick_params(labelsize=cbar_ticklab_fontsize)
            cax.set_ylabel(cbar_axlab, fontdict={'fontsize': cbar_axlab_fontsize})
            cax.yaxis.get_offset_text().set_fontsize(cbar_ticklab_fontsize)
            cax.yaxis.set_offset_position('right')
        # format axes
        ax.set_xlim(rast.rio.bounds()[0::2])
        ax.set_ylim(rast.rio.bounds()[1::2])
        phf.strip_axes_labels_and_ticks(ax)
        if row == 0:
            ax.set_title('%i km neighborhood' % neigh_rad,
                         fontdict={'fontsize': 9})
        del rast



fig_supp = plt.figure(figsize=(8.75,10))
gs = fig_supp.add_gridspec(nrows=7,
                           ncols=5,
                           wspace=0,
                           hspace=0,
                           width_ratios=[1,1,1,0.08, 0.02],
                          )

# make vars' supp figs (each one stacking all 3 neighborhood radii)
for n, var in enumerate(vars):
    print('\n\nNOW PLOTTING %s..\n\n' % var)
    map_asynch(fig=fig_supp,
               gs=gs,
               row=n,
               cbar_axlab=cbar_axlab_dict[var],
               var=var,
               cbar_axlab_fontsize=cbar_axlab_fontsize,
               cbar_ticklab_fontsize=cbar_ticklab_fontsize,
              )
fig_supp.subplots_adjust(bottom=0.01, top=0.97, left=0.01, right=0.88, hspace=0, wspace=0)
fig_supp.savefig(os.path.join(phf.FIGS_DIR, 'FIG_SUPP_all_asynch_maps.png'), dpi=600)



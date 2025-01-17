import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
from rasterstats import zonal_stats
from rasterio.enums import Resampling
import os
import sys
import re


sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf

# set data-directory path
data_dir = phf.EXTERNAL_DATA_DIR

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



def map_r2(ax, r2_filename, axlabel,
           add_colorbar=False,
           cbar_ax=None,
           cbar_kwargs=None,
          ):
    """
    plot an R2 map from the LSP/seasonality-fitting process
    """
    files = [f for f in os.listdir(data_dir) if
             os.path.split(f)[-1] == r2_filename]
    assert len(files) == 1
    r2 = rxr.open_rasterio(os.path.join(data_dir, files[0]), masked=True)[0]
    # NOTE: annoying AttributeError is because da.attrs['long_name']
    #       is retained as a tuple of names (rather than being subsetted
    #       by indexing) when I index a single layer out of an
    #       xarray.core.dataarray.DataArray;
    #       for now, a hacky fix is just assigning a string to that attr
    r2.attrs['long_name'] = ''
    cmap = 'viridis'
    r2.plot.imshow(ax=ax,
                   cmap=cmap,
                   vmin=0,
                   vmax=1,
                   add_colorbar=add_colorbar,
                   cbar_ax=cbar_ax,
                   cbar_kwargs=cbar_kwargs,
                   zorder=0,
                  )
    phf.plot_juris_bounds(ax,
                          lev0_linecolor='gray',
                          lev0_linewidth=0.5,
                          lev0_alpha=0.8,
                          lev0_zorder=2,
                          lev1_linecolor='gray',
                          lev1_linewidth=0.3,
                          lev1_alpha=0.7,
                          lev1_zorder=1,
                          crs=r2.rio.crs.to_epsg(),
                          strip_axes=True,
                         )
    ax.text(-165, -75, axlabel, size=9, clip_on=False)


if __name__ == '__main__':

    # dict for labeling maps
    label_dict = {'NIRv': 'NIRv',
                  'SIF': 'SIF',
                  'pr': 'mean ann. precip.',
                  'tmmn': 'mean ann. min. monthly temp.',
                  'tmmx': 'mean ann. max. monthly temp.',
                  'def': 'mean ann. climate water deficit',
                  'cloud': 'cloud-cover percent',
                 }

    fig = plt.figure(figsize=(8.7, 8))
    gs = GridSpec(4, 2, figure=fig)
    filenames = ['NIRv_harm_reg_R2.tif',
                 'SIF_harm_reg_R2.tif',
                 'pr_harm_reg_R2.tif',
                 'tmmn_harm_reg_R2.tif',
                 'tmmx_harm_reg_R2.tif',
                 'def_harm_reg_R2.tif',
                 'cloud_harm_reg_R2.tif',
                ]
    for ct, fn in enumerate(filenames):
        print(f"\n\nNOW PROCESSING {fn}...\n\n")
        ax = fig.add_subplot(gs[ct//2, ct%2])
        label = label_dict[fn.replace('_harm_reg_R2.tif', '')]
        if ct == len(filenames)-1:
            add_colorbar = True
            cbar_ax=fig.add_subplot(GridSpec(100, 100, figure=fig)[88:92,60:90])
            cbar_kwargs = {'orientation': 'horizontal'}
        else:
            add_colorbar = False
            cbar_ax = None
            cbar_kwargs = None
        map_r2(ax,
               fn,
               label,
               add_colorbar=add_colorbar,
               cbar_ax=cbar_ax,
               cbar_kwargs=cbar_kwargs,
              )
    fig.axes[-1].set_xlabel('$R^2$', fontdict={'fontsize': 10})
    fig.subplots_adjust(left=0.05,
                        right=0.95,
                        bottom=0.01,
                        top=0.99,
                        hspace=0,
                        wspace=0.2,
                       )
    fig.savefig(os.path.join(phf.FIGS_DIR, 'FIG_SUPP_harm_reg_R2_maps.png'), dpi=500)


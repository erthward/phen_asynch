import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
from shapely.geometry import Point
from shapely.geometry import Polygon as shapelyPolygon
import numpy as np
import xarray as xr
import rasterio as rio
import rioxarray as rxr
import re, os, sys

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                                                            'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf


# data directories
rs_datadir = phf.EXTERNAL_DATA_DIR
flux_datadir = phf.EXTERNAL_FLUX_DATA_DIR

# load NIRv-SIF R2s map
r2s = rxr.open_rasterio(os.path.join(rs_datadir, 'NIRv_SIF_phen_R2s.tif'),
                        masked=True)[0]
print(f'\n\nMEDIAN NIRv-SIF R2 VALUE: {np.nanmedian(r2s)}\n\n')

# map results
fig = plt.figure(figsize=(7,4.6))
ax = fig.add_subplot(1, 1, 1)
divider = make_axes_locatable(ax)
cax = divider.append_axes('bottom', size='8%', pad=0.5)
r2s.plot.imshow(ax=ax,
                vmin=0,
                vmax=1,
                cmap='gray',
                zorder=2,
                add_colorbar=True,
                cbar_ax=cax,
                cbar_kwargs={'orientation': 'horizontal'},
               )
phf.plot_juris_bounds(ax,
                      lev0_color='#ede6d1', #'#9e8e67',
                      lev0_alpha=0.8,
                      lev0_zorder=0,
                      lev1_linewidth=0.25,
                      lev1_alpha=0.6,
                      lev1_zorder=1,
                      crs = r2s.rio.crs.to_epsg(),
                      strip_axes=True,
                     )
ax.set_xlim(r2s.rio.bounds()[::2])
ax.set_ylim(r2s.rio.bounds()[1::2])
ax.text(1.14 * r2s.rio.bounds()[0],
        0.96 * r2s.rio.bounds()[3],
        'A.',
        weight='bold',
        size=20,
       )
cax.set_xlabel('')
cax.text(0.49, 1.25, '$R^2$', size=16)

fig.subplots_adjust(top=0.98, bottom=0.04, left=0.05, right=0.97)

fig.savefig(os.path.join(phf.FIGS_DIR, 'FIG_SUPP_NIRv-SIF_LSP_comparison_results.png'), dpi=700)


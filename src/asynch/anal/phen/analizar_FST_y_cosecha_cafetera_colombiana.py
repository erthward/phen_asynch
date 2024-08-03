#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
import palettable
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
from sklearn.cluster import KMeans
from shapely.geometry import Point
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf
from MMRR import MMRR


# load digitized cosecha points
patt = '(VERDE)|(AMARILLO)|(ANARANJADO)|(MORADO)'
files = [f for f in os.listdir('.') if re.search(patt, f)]
# reorder files so that they plot in the overal N-S order of the clusters
NS_color_order = ['VERDE', 'ANARANJADO', 'MORADO', 'AMARILLO']
order = [np.argwhere([re.search(c, f) for f in files]) for c in NS_color_order]
order = np.stack(order).ravel()
files = [*np.array(files)[order]]
# I've named files based on the colors used in the figure I digitized, but we
# should instead choose colors that better highlight the overarching N-S
# shift in LSP shape and timing (but that are still distinct enough to clearly
# highlight the very short-distance regional discontinuities)
cmap = mpl.cm.viridis
rgbas = cmap(np.linspace(0, 1, 4))
hexes = [mpl.colors.rgb2hex(v, keep_alpha=True) for v in rgbas]
colors = {'VERDE': hexes[3],
          'ANARANJADO': hexes[2],
          'MORADO': hexes[1],
          'AMARILLO': hexes[0],
         }
gdfs = []
for f in files:
    df = pd.read_csv(f, header=None)
    df.loc[:, 'geometry'] = [Point(*row.values) for i, row in df.iterrows()]
    principal = re.search('(?<=principal)[A-Z]{4}', f).group()
    mitaca = re.search('(?<=mitaca)[A-Z]{2}', f)
    df.loc[:, 'principal'] = principal
    if mitaca is not None:
        df.loc[:, 'mitaca'] = mitaca.group()
    else:
        df.loc[:, 'mitaca'] = np.nan
    color = re.search(patt, f).group()
    df.loc[:, 'color'] = colors[color]
    gdf = gpd.GeoDataFrame(df, geometry='geometry', crs=4326)
    gdfs.append(gdf)
gdf = pd.concat(gdfs)
# reset index, to facilitate dropping points with missing LSP data
gdf = gdf.reset_index()

# set up figure
fig = plt.figure(figsize=(7, 5))
gs = fig.add_gridspec(nrows=50, ncols=10)

# NOTE: numerical doy range associated with the major and minor seasons,
#       adjusted down by 1 to account for January 1st being doy 0
season_xlims = {'MAMJ': [59, 180],
                'SOND': [243, 364],
                'AM': [90, 150],
                'ON': [273, 333],
               }
all_keep = []
all_lsp = []
for i, color in enumerate(colors.values()):
    ax = fig.add_subplot(gs[11*i:11*(i+1), 5:])
    # extract LSP ts at each point
    subgdf = gdf.loc[gdf['color'] == color]
    pts = subgdf.get_coordinates().values
    lsp_ts = phf.get_raster_info_points(phf.COEFFS_FILE,
                                        pts,
                                        'ts',
                                        standardize=True,
                                        fill_nans=False,
                                        fill_tol=None,
                                       )
    assert np.all(lsp_ts.shape == (len(pts), 365))
    # drop NA points
    keep = np.sum(pd.isnull(lsp_ts), axis=1) == 0
    all_keep.extend([*subgdf.index.values[keep]])
    print((f"\n\t{np.round(100*np.mean(np.invert(keep)), 1)}% of points "
           f"in cluster {i} are missing LSP data...\n"))
    lsp_ts = lsp_ts[keep, :]
    all_lsp.append(lsp_ts)
    mid = np.median(lsp_ts, axis=0)
    # plot regions' low and high percentiles
    pctile_lo = np.percentile(lsp_ts, 10, axis=0)
    pctile_hi = np.percentile(lsp_ts, 90, axis=0)
    pctile_poly_coords = np.vstack((np.array([*zip(range(365), pctile_lo)]),
                        np.array([*zip([*range(365)][::-1], pctile_hi[::-1])])))
    pctile_poly = Polygon(pctile_poly_coords)
    patchcoll = PatchCollection([pctile_poly],
                                alpha=0.3,
                                color=color,
                                edgecolor='k',
                               )
    ax.add_collection(patchcoll)
    # plot the median as a darker line on top
    ax.plot(range(len(mid)),
            mid,
            linestyle='-',
            color=color,
            alpha=1,
            zorder=1,
            linewidth=3,
           )
    principal = np.unique(subgdf['principal'])
    assert len(principal) == 1
    principal = principal[0]
    ax.plot(season_xlims[principal],
            [-2.3]*2,
            color=color,
            linewidth=5,
            zorder=0,
           )
    if np.sum(pd.isnull(subgdf['mitaca'])) == 0:
        mitaca = np.unique(subgdf['mitaca'])
        assert len(mitaca) == 1
        mitaca = mitaca[0]
        ax.plot(season_xlims[mitaca],
                [-2.3]*2,
                color=color,
                linewidth=5,
                zorder=0,
               )

    ax.set_xlim(0, 364)
    ax.set_ylim(-2.5, 2)
    ax.set_yticks(())
    if i == 2:
        ax.set_ylabel(f"{' '*18}scaled LSP", fontdict={'fontsize': 14})
    else:
        ax.set_ylabel('')
    if i == 3:
        ax.set_xlabel('day of year', fontdict={'fontsize': 14})
        ax.set_xticks([0, 90, 181, 273, 364],
                      ['Jan', 'Apr', 'Jul', 'Oct', 'Jan'],
                     )
        ax.tick_params(labelsize=10)
    else:
        ax.set_xlabel('')
        ax.set_xticks(())

# keep only points with LSP data
gdf = gdf.loc[all_keep, :]

# read in the unfolded EOFs file
# (NOTE: avoids the 'color smuding' artefact caused by folding the EOFS across
# the ITCZ, which isn't an issue in most other focal/regional maps, but
# is a concern here)
eofs = rxr.open_rasterio(phf.EOFS_FILE)[:3].rio.set_crs(4326)
for i in range(3):
    eofs[i] = (eofs[i] - np.nanmin(eofs[i]))/(np.nanmax(eofs[i])-np.nanmin(eofs[i]))
eofs = eofs.rio.reproject(8857)
gdf_plot = gdf.to_crs(8857)

# map regions, as concave hulls, on top of EOFs map
ax = fig.add_subplot(gs[:49, :5])
ax.set_aspect('equal')
eofs.plot.imshow(ax=ax)
for color in np.unique(gdf['color']):
    subgdf = gdf_plot[gdf_plot['color'] == color]
    subgdf.plot(ax=ax,
                color=color,
                edgecolor='k',
                marker='.',
                markersize=28,
                linewidth=0.1,
                zorder=2,
               )

phf.plot_juris_bounds(ax,
                      lev0_linewidth=0.3,
                      lev0_zorder=1,
                      lev1=False,
                      crs=8857,
                      strip_axes=True,
                      reset_axlims=False,
                     )
ax.set_xlim(-7.45e6, -6.863e6)
ax.set_ylim(1.38e5, 1.46e6)

# format and save figure
fig.subplots_adjust(hspace=0,
                    left=0,
                    right=0.98,
                    bottom=0,
                    top=0.98,
                   )
fig.savefig(os.path.join(phf.FIGS_DIR,
                         'FST_y_cosecha_cafetero_colombiana.png'),
            dpi=600,
           )


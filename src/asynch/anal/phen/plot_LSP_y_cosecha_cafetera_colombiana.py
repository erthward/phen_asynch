#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from shapely.geometry import Point
from matplotlib.colors import ListedColormap

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf
from MMRR import MMRR



# load digitized cosecha points
patt = '(VERDE)|(AMARILLO)|(ANARANJADO)|(MORADO)'
files = [f for f in os.listdir('.') if re.search(patt, f)]
colors = {'AMARILLO': 'yellow',
          'VERDE': 'limegreen',
          'ANARANJADO': 'orange',
          'MORADO': 'purple',
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


# plot regions' median or mean LSP curves
fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(2, 1, 2)
mid_fn = np.median
# NOTE: numerical doy range associated with the major and minor seasons,
#       adjusted down by 1 to account for January 1st being doy 0
season_xlims = {'MAMJ': [59, 180],
                'SOND': [243, 364],
                'AM': [90, 150],
                'ON': [273, 333],
               }
all_keep = []
all_lsp = []
for i, color in enumerate(np.unique(gdf['color'])):
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
    # drop NA points
    keep = np.sum(pd.isnull(lsp_ts), axis=1) == 0
    all_keep.extend([*subgdf.index.values[keep]])
    print((f"\n\t{np.round(100*np.mean(np.invert(keep)), 1)}% of points "
           f"in {color} cluster are missing LSP data...\n"))
    lsp_ts = lsp_ts[keep, :]
    all_lsp.append(lsp_ts)
    mid = mid_fn(lsp_ts, axis=0)
    pctile_05 = np.percentile(lsp_ts, 5, axis=0)
    pctile_95 = np.percentile(lsp_ts, 95, axis=0)
    ax.plot(range(len(mid)),
            mid,
            linestyle='-',
            color='k',
            alpha=1,
            zorder=0,
            linewidth=0.5,
           )
    ax.plot(range(len(mid)),
            mid,
            linestyle='-',
            color=color,
            alpha=0.8,
            zorder=1,
            linewidth=3,
           )
    ax.plot(range(len(mid)),
            pctile_05,
            linestyle='--',
            linewidth=0.5,
            color=color,
            alpha=0.7,
            zorder=1,
           )
    ax.plot(range(len(mid)),
            pctile_95,
            linestyle='--',
            linewidth=0.5,
            color=color,
            alpha=0.7,
            zorder=1,
           )
    principal = np.unique(subgdf['principal'])
    assert len(principal) == 1
    principal = principal[0]
    ax.plot(season_xlims[principal],
            [-2-(i*0.1)]*2,
            color=color,
            linewidth=5,
            zorder=0,
           )
    if np.sum(pd.isnull(subgdf['mitaca'])) == 0:
        mitaca = np.unique(subgdf['mitaca'])
        assert len(mitaca) == 1
        mitaca = mitaca[0]
        ax.plot(season_xlims[mitaca],
                [-2-(i*0.1)]*2,
                color=color,
                linewidth=5,
                zorder=0,
               )
ax.set_xlim(0, 364)
ax.set_ylim(-2.5, 2)
ax.set_ylabel('scaled LSP', fontdict={'fontsize': 16})
ax.set_yticks(())
ax.set_xlabel('day of year', fontdict={'fontsize': 16})
ax.set_xticks([0, 90, 181, 273, 364],
              ['Jan', 'Apr', 'Jul', 'Oct', 'Jan'],
             )
ax.tick_params(labelsize=12)

# keep only points with LSP data
gdf = gdf.loc[all_keep, :]

# read in prepped EOFs file
eofs = rxr.open_rasterio(phf.EOFS_PREPPED_FILE)[:3].rio.set_crs(8857)
eofs = eofs.rio.reproject(4326)

# map regions, as concave hulls, on top of EOFs map
ax = fig.add_subplot(2, 2, 1)
eofs.plot.imshow(ax=ax)
#gdf.plot(ax=ax,
#         color='white',
#         edgecolor='black',
#         markersize=5,
#         zorder=1,
#        )
for color in np.unique(gdf['color']):
    subgdf = gdf[gdf['color'] == color]
    subgdf.dissolve().concave_hull().plot(ax=ax,
                                          color=color,
                                          linewidth=0.1,
                                          zorder=1,
                                         )
    #subgdf.plot(ax=ax,
    #            color=subgdf['color'],
    #            markersize=8,
    #            marker='.',
    #            zorder=1,
    #           )

phf.plot_juris_bounds(ax,
                      lev1_alpha=0.6,
                      lev1_linewidth=0.05,
                      lev1_zorder=2,
                      lev0_linewidth=0.1,
                      lev0_zorder=3,
                      crs=4326,
                      strip_axes=True,
                      reset_axlims=False,
                     )
ax.set_xlim(-79, -70)
ax.set_ylim(1, 12)


# map locations, colored by LSP clustering
ax = fig.add_subplot(2, 2, 2)
lsp_gdf = gdf.loc[all_keep, ['geometry', 'color']]
all_lsp = np.vstack(all_lsp)
for i in range(365):
    lsp_gdf.loc[:, f'doy{i}'] = all_lsp[:, i]
# scree plot suggests k=3 is best
np.random.seed(1)
lsp_clusts = KMeans(n_clusters=3).fit(lsp_gdf.iloc[:, 2:].values)
lsp_gdf.loc[:, 'clust'] = lsp_clusts.labels_
lsp_gdf.plot(ax=ax,
             column='clust',
             markersize=8,
             marker='.',
             cmap=ListedColormap(['limegreen', 'yellow', 'orange']),
            )
phf.plot_juris_bounds(ax,
                      lev1_alpha=0.6,
                      lev1_linewidth=0.05,
                      lev1_zorder=1,
                      lev0_linewidth=0.1,
                      lev0_zorder=2,
                      crs=4326,
                      strip_axes=True,
                      reset_axlims=False,
                     )
ax.set_xlim(-79, -70)
ax.set_ylim(1, 12)

fig.savefig('LSP_y_cosecha_cafetero_colombiana.png',
            dpi=600,
           )


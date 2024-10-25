#!/usr/bin/env python
# coding: utf-8

import os
import re
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
import matplotlib.pyplot as plt
from shapely.geometry import Point
from datetime import datetime, timedelta

import sys
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf


# load the extracted data
extract = gpd.read_file(os.path.join(phf.EXTERNAL_DATA_DIR,
                                     'raw_NIRv_at_CA_sites.geojson'))

# clean it up
date_strs = [re.search('\d{4}_\d{2}_\d{2}', s).group() for s in extract['id']]
extract.loc[:, 'date'] = [datetime.strptime(s, '%Y_%m_%d') for s in date_strs]
NIRv_dfs = {}
locs = {'site': [],
        'geometry': [],
       }
for site in extract['SITE_NAME'].unique():
    df = extract[extract['SITE_NAME'] == site]
    unique_locs = df['geometry'].unique()
    assert len(unique_locs) == 1
    geom = df['geometry'].iloc[0]
    coords = (geom.xy[0][0], geom.xy[1][0])
    locs['site'].append(site)
    locs['geometry'].append(Point(coords))
    df = df.loc[:, ['date', 'mean']]
    df = df.rename(mapper = {'mean': 'NIRv'}, axis=1)
    df.sort_values(by='date', inplace=True)
    NIRv_dfs[site] = df

# make locations GeoDataFrame
locs_df = pd.DataFrame(locs)
locs_gdf = gpd.GeoDataFrame(locs_df, crs=4326, geometry='geometry')

# common equal-area projection to use
crs = 4326

# prep eofs file to use for plot
eofs = rxr.open_rasterio(phf.EOFS_PREPPED_FILE)[:3]
eofs.rio.set_crs(8857)
eofs = eofs.rio.reproject(crs)

# set the map limits
xlims = (-123.5, -119)
ylims = (36.5, 40.5)

# set up the figure
fig = plt.figure(figsize=(23,10))
gs = fig.add_gridspec(nrows=100, ncols=230)
ax_rgb = fig.add_subplot(gs[:, 3:80])
ts_axs = [fig.add_subplot(gs[i0:i1, 90:]) for i0, i1 in zip([2, 33, 64],
                                                            [33, 64, 95])]

# plot the EOFS map and the sites
eofs.plot.imshow(ax=ax_rgb)
ax_rgb.set_xlim(xlims)
ax_rgb.set_ylim(ylims)
phf.plot_juris_bounds(ax_rgb,
                      crs=crs,
                      strip_axes=False,
                      reset_axlims=False,
                     )
locs_gdf.plot(color='black',
              markersize=30,
              ax=ax_rgb,
             )
for i, row in locs_gdf.iterrows():
    ax_rgb.text(row.geometry.x+0.023,
                row.geometry.y+0.01,
                row['site'],
                rotation=0,
                size=18,
                weight='bold',
               )
# plot Cal, as an easter egg! :)
ax_rgb.scatter(-122.25779,
               37.87235,
               marker='.',
               color='#271a69',
               s=2,
              )
ax_rgb.set_xlabel('longitude', fontdict={'size': 17})
ax_rgb.set_ylabel('latitude', fontdict={'size': 17})
ax_rgb.set_xticks(ax_rgb.get_xticks(), ax_rgb.get_xticklabels(), rotation=45)
ax_rgb.tick_params(labelsize=11)
ax_rgb.set_title('')

# get and plot both the raw and fitted time series for each site
# NOTE: divide values by 10000, since I didn't do that in GEE,
#       at some point in an effort to use ints to reduce memory usage
i = 0
for site, df in NIRv_dfs.items():
    # get the fitted LSP ts
    pts = np.stack(locs_df.loc[locs_df['site'] == site, 'geometry'].values[0].xy).T
    NIRv_fit = phf.get_raster_info_points(phf.COEFFS_FILE,
                                          pts,
                                          'ts',
                                          standardize=False,
                                          fill_nans=False,
                                          fill_tol=None,
                                         ).ravel()
    # reconcile fitted to raw dates
    doys = [d.timetuple().tm_yday for d in df['date']]
    # NOTE: just clip last day of leap year to 365, for simplicity sake
    doys = np.clip(doys, a_min=1, a_max=365)
    NIRv_fit_rec = np.array([NIRv_fit[doy-1] for doy in doys])
    # plot raw and fitted
    ax_ts = ts_axs[i]
    ax_ts.plot(df['date'],
               df['NIRv']/10000,
               '-k',
               zorder=1,
               label='raw',
              )
    ax_ts.plot(df['date'],
               NIRv_fit_rec/10000,
               '--r',
               zorder=2,
               label='fitted',
              )
    if i == 0:
        ax_ts.legend(fontsize=16)
    ax_ts.set_title(site,
                    y=0.88,
                    fontdict={'size': 19,
                              'weight': 'bold',
                             },
                   )
    ax_ts.set_ylabel('$NIR_{V}$',
                     fontdict={'fontsize': 17},
                    )
    # add vertical gray lines at each year
    years = pd.date_range(df['date'].iloc[0],
                          df['date'].iloc[-1]+timedelta(days=1),
                          21,
                         )
    for yr in years:
        ax_ts.axvline(x=yr,
                      ymin=0,
                      ymax=1,
                      linestyle=':',
                      color='gray',
                      alpha=0.5,
                      zorder=0,
                     )
    ax_ts.set_xlim(years[0], years[-1])
    # NOTE: set ylim to a value comfortably above the max across the datasets,
    #       rather than 1, to make sure temporal variability is clearly visible
    ax_ts.set_ylim(0, 0.3)
    if i > 0:
        ax_ts.set_yticks([0, 0.1, 0.2])
    else:
        ax_ts.set_yticks([0, 0.1, 0.2, 0.3])
    if i == 2:
        ax_ts.set_xticks(years)
        ax_ts.set_xticklabels([str(yr.year) for yr in years],
                              rotation=45,
                             )
        ax_ts.set_xlabel('date',
                         fontdict={'fontsize': 17},
                        )
    else:
        ax_ts.set_xticks(())
        ax_ts.set_xlabel('')
    ax_ts.tick_params(labelsize=14)
    i+=1

# adjust subplots and save
fig.subplots_adjust(hspace=0,
                    wspace=0,
                    left=0.03,
                    right=0.98,
                    bottom=0.05,
                    top=1.0,
                   )
fig.savefig(os.path.join(phf.FIGS_DIR, 'FIG_SUPP_raw_and_fitted_NIRv_LSP.png'),
            dpi=600,
           )


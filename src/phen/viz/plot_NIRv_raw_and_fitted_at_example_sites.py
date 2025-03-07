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
from shapely.geometry import Point
from datetime import datetime, timedelta

import sys
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf

# get the site pair name to plot
try:
    spn = sys.argv[1]
except Exception:
    spn = input('spn: ')
assert spn in ['boreal', 'andes', 'outback']

# load the extracted data
extract = gpd.read_file(os.path.join(phf.EXTERNAL_DATA_DIR,
                                     'raw_NIRv_at_test_sites.geojson'))

# clean it up
date_strs = [re.search('\d{4}_\d{2}_\d{2}', s).group() for s in extract['id']]
extract.loc[:, 'date'] = [datetime.strptime(s, '%Y_%m_%d') for s in date_strs]
NIRv_dfs = {}
locs = {'site': [],
        'geometry': [],
       }
for site in extract['name'].unique():
    df = extract[extract['name'] == site]
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

# strings identifying site pairs
site_pair_names = ['boreal', 'andes', 'outback']

# make locations GeoDataFrame
locs_df = pd.DataFrame(locs)
locs_gdf = gpd.GeoDataFrame(locs_df, crs=4326, geometry='geometry')

# common equal-area projection to use
crs = 4326

# map limits for each pair of sites
xlims = {'boreal': (-104.05, -102.85),
         'andes': (-74.49, -73.11),
         'outback': (129.05, 132.45),
        }
ylims = {'boreal': (52.70, 53.27),
         'andes': (4.96 , 5.88),
         'outback': (-32.02, -27.88),
        }

# prep eofs file to use for plot
eofs = rxr.open_rasterio(phf.EOFS_PREPPED_FILE)[:3]
eofs.rio.set_crs(8857)
eofs = eofs.rio.reproject(crs)

# set up the figure
fig = plt.figure(figsize=(23,8))
gs = fig.add_gridspec(nrows=100, ncols=230)
ax_rgb = fig.add_subplot(gs[:, 3:120])
ts_axs = [fig.add_subplot(gs[i0:i1, 130:]) for i0, i1 in zip([2, 52],
                                                            [50, 100])]

# plot the EOFS map and the sites
eofs.plot.imshow(ax=ax_rgb)
ax_rgb.set_xlim(xlims[spn])
ax_rgb.set_ylim(ylims[spn])
phf.plot_juris_bounds(ax_rgb,
                      crs=crs,
                      strip_axes=False,
                      reset_axlims=False,
                     )
locs_subgdf = locs_gdf[[site.startswith(spn) for site in locs_gdf['site']]]
locs_subgdf.plot(color='black',
                 markersize=30,
                 ax=ax_rgb,
                )
for i, row in locs_subgdf.iterrows():
    ax_rgb.text(row.geometry.x+0.023,
                row.geometry.y+0.01,
                row['site'].replace('_', ' ').capitalize(),
                rotation=0,
                size=18,
                weight='bold',
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
    if site.startswith(spn):
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
        ax_ts.set_title(site.replace('_', ' ').capitalize(),
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
        max_max = max([np.nanmax(df['NIRv'].values) for k,
                            df in NIRv_dfs.items() if k.startswith(spn)])
        ax_ts.set_ylim(0, (max_max/10000)+0.1)
        ax_ts.set_yticks(np.arange(0, max_max/10000 + 0.1, .1))
        if i == 1:
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
                    bottom=0.12,
                    top=0.97
                   )
fig.savefig(os.path.join(phf.FIGS_DIR, f'FIG_SUPP_raw_and_fitted_NIRv_LSP_{spn}.png'),
            dpi=600,
           )

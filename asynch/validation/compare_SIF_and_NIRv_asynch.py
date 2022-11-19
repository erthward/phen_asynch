import numpy as np
import rioxarray as rxr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import os, sys, re

# data directory
data_dir = "/media/deth/SLAB/diss/3-phn/GEE_outputs/final/"
countries_data_dir = "/home/deth/Desktop/CAL/research/projects/seasonality/seasonal_asynchrony/data/bounds/"

# load country boundaries
countries = gpd.read_file(os.path.join(countries_data_dir, 'NewWorldFile_2020.shp'))


# loop over neighborhood radii (in km)
fig = plt.figure(figsize=(34,30))
gs = fig.add_gridspec(nrows=3, ncols=4, width_ratios=[1,1,0.1,1])
for neigh_rad_i, neigh_rad in enumerate(['50', '100', '150']):

    print('\n\nRUNNING COMPARISON FOR %s KM-RADIUS NEIGHBORHOOD...\n\n' % neigh_rad)

    # load both rasters
    sif = rxr.open_rasterio(os.path.join(data_dir,
                'SIF_STRICT_asynch_%skm.tif' % neigh_rad), masked=True)[0]
    nirv = rxr.open_rasterio(os.path.join(data_dir,
                'NIRv_STRICT_asynch_%skm.tif' % neigh_rad), masked=True)[0]

    # center and scale each one ((x-mu)/sigma)
    sif_scale = (sif-np.nanmean(sif))/np.nanstd(sif)
    nirv_scale = (nirv-np.nanmean(nirv))/np.nanstd(nirv)

    # get difference and plot
    diff_scale = nirv_scale - sif_scale

    # get max abs val to use to center raster map on 0
    max_abs_val = max(np.abs((np.nanpercentile(diff_scale, 0.01),
                              np.nanpercentile(diff_scale, 99.99))))

    ax = fig.add_subplot(gs[neigh_rad_i, :2])
    countries.to_crs(diff_scale.rio.crs).plot(ax=ax,
                                              color='none',
                                              edgecolor='black',
                                              zorder=1)
    diff_scale.plot.imshow(ax=ax,
                           cmap='coolwarm_r',
                           vmin=-max_abs_val,
                           vmax=max_abs_val,
                           zorder=0)
    # increase colorbar ticklabel size and label the colorbar
    fig.axes[-1].tick_params(labelsize=20)
    fig.axes[-1].set_ylabel('$NIR_{V\ stand}-SIF_{stand}$',
                            fontdict={'fontsize': 30})

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_title('neigh_rad = %s km' % neigh_rad,
                 fontdict={'fontsize': 45})

    # scatter samples against one another and fit SLR
    ax = fig.add_subplot(gs[neigh_rad_i, 3])
    ax.scatter(sif.values.ravel(),
               nirv.values.ravel(),
               s=0.1,
               c='black',
               alpha=0.01,
              )

    endog = np.atleast_1d(nirv.values.ravel())
    exog = np.atleast_1d(sif.values.ravel())
    endog = endog[pd.notnull(exog)]
    exog = exog[pd.notnull(exog)]
    exog = exog[pd.notnull(endog)]
    endog = endog[pd.notnull(endog)]
    lm = sm.OLS(endog=endog.T, exog=np.stack((np.ones(exog.shape), exog)).T).fit()

    plot_x = np.linspace(np.nanpercentile(sif.values, 0.01),
                         np.nanpercentile(sif.values, 99.99),
                         10000)
    plot_y = lm.params[0] + lm.params[1]*plot_x
    ax.plot(plot_x, plot_y, '-r', linewidth=2, alpha=0.5)
    ax.set_xlim(np.nanpercentile(sif, (1, 99)))
    ax.set_ylim(np.nanpercentile(nirv, (1, 99)))
    get_text_pos = lambda lims, frac: lims[0] + (frac*(lims[1]-lims[0]))
    ax.text(get_text_pos(ax.get_xlim(), 0.8),
            get_text_pos(ax.get_ylim(), 0.1),
            '$R^{2}=%0.2f$' % lm.rsquared,
            color='red',
            size=26,
           )
    ax.set_xlabel('$SIF$ asynchrony', fontdict={'fontsize': 30})
    ax.set_ylabel('$NIR_{V}$ asynchrony', fontdict={'fontsize': 30})
    ax.tick_params(labelsize=20, rotation=45)

fig.subplots_adjust(top=0.95,
                    bottom=0.08,
                    left=0.02,
                    right=0.96,
                    hspace=0.2,
                    wspace=0)
fig.savefig('FIG_S4_scaled_NIRv_vs_scaled_SIF.png', dpi=700)

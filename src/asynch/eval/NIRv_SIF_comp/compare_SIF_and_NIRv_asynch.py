import numpy as np
import rioxarray as rxr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import os, sys, re

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                                                            'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf


# data directory
data_dir = phf.EXTERNAL_DATA_DIR

# loop over neighborhood radii (in km)
fig = plt.figure(figsize=(12,10))
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
    max_abs_val = max(np.abs((np.nanpercentile(diff_scale, 1),
                              np.nanpercentile(diff_scale, 99))))

    ax = fig.add_subplot(gs[neigh_rad_i, :2])
    diff_scale.plot.imshow(ax=ax,
                           cmap='coolwarm_r',
                           vmin=-max_abs_val,
                           vmax=max_abs_val,
                           zorder=0)
    # increase colorbar ticklabel size and label the colorbar
    fig.axes[-1].tick_params(labelsize=9)
    fig.axes[-1].set_ylabel('$asynch_{NIR_{V\ stand}}-asynch_{SIF_{stand}}$',
                            fontdict={'fontsize': 14})

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_title('neighborhood radius = %s km' % neigh_rad,
                 fontdict={'fontsize': 16})
    phf.plot_juris_bounds(ax=ax,
                          crs=diff_scale.rio.crs.to_epsg(),
                          strip_axes=False,
                         )
    # truncate to match northern extent of NIRv dataset
    phf.set_upper_ylim(ax, uplim=60)

    # scatter samples against one another and fit SLR
    ax = fig.add_subplot(gs[neigh_rad_i, 3])
    ax.ticklabel_format(style='sci', scilimits=([-5,-5]))
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
    ax.set_xlim(np.nanpercentile(sif, (0.1, 99.9)))
    ax.set_ylim(np.nanpercentile(nirv, (0.1, 99.9)))
    get_text_pos = lambda lims, frac: lims[0] + (frac*(lims[1]-lims[0]))
    ax.text(get_text_pos(ax.get_xlim(), 0.7),
            get_text_pos(ax.get_ylim(), 0.8),
            '$R^{2}=%0.2f$' % lm.rsquared,
            color='red',
            size=11,
           )
    if neigh_rad_i == 2:
        ax.set_xlabel('$SIF$ asynchrony', fontdict={'fontsize': 14})
    else:
        ax.set_xlabel('')
    ax.set_ylabel('$NIR_{V}$ asynchrony', fontdict={'fontsize': 14})
    ax.tick_params(labelsize=11, rotation=0)

fig.subplots_adjust(top=0.95,
                    bottom=0.08,
                    left=0.02,
                    right=0.96,
                    hspace=0.2,
                    wspace=0.2,
                   )
fig.savefig(os.path.join(phf.FIGS_DIR,
                         'FIG_SUPP_scaled_NIRv_asynch_vs_scaled_SIF_asynch.png'), dpi=500)


import os
import re
import sys
import numpy as np
import pandas as pd
import rasterio as rio
import geopandas as gpd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from shapely.geometry import Point


#------------------------------------------------------------
# TODO:

# what to do about very late 'first's?
'''
fig, ax = plt.subplots(1,1)
nam.plot(color='none', ax=ax)
npn[npn['mean_first_yes_doy']>300].plot('mean_first_yes_doy', ax=ax, zorder=1)
print(npn[npn['mean_first_yes_doy']>300].loc[:, ['genus', 'species']])
'''

# drop duplicate sites and deal with spatial autocorrelation!

# should SI-x regression be logistic instead?

#------------------------------------------------------------




# local modules
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf

# whether to interpolate missing LSP coefficients, 
# and if so, how far away can rasterio look for neighbor to interpolate from?
# coefficient raster?
interp_sea_data = False
neigh_dist_sea_fill_tol = 5

# load NPN data
npn = pd.read_csv('./NPN_leaves_data_dom_tree_spp.csv')
npn['geometry'] = [Point(row['longitude'],
                         row['latitude']) for i, row in npn.iterrows()]
gpd.GeoDataFrame(npn, geometry='geometry', crs=4326)
npn = gpd.GeoDataFrame(npn, geometry='geometry', crs=4326)

# load basic NAm shapes
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
nam = world[world['continent'] == 'North America']
del world

# NOTE: needs to be based on standardized time series,
#       to focus on timing rather than fitted NIRV values
pts = np.array([[g.x, g.y] for g in npn['geometry'].values])
lsp = phf.get_raster_info_points(phf.COEFFS_FILE,
                                 pts,
                                 'ts',
                                 standardize=True,
                                 fill_nans=interp_sea_data,
                                 fill_tol=neigh_dist_sea_fill_tol,
                                )
print((f"\n\n{np.round(100*np.mean(np.any(pd.isnull(lsp), axis=1)), 2)} % "
        "of sea_dist consists of missing values"))

# subset to points where we have data
not_missing = np.all(pd.notnull(sea_dist), axis=1)
npn = npn.loc[not_missing, ]
lsp = lsp[not_missing, :]
assert lsp.shape[0] == len(npn)

# estimate start of season (SOS) as doy when LSP exceeds 50% of its amplitude
# (a la Melaas et al. 2016, but I take it this is standard practice)
half_max = np.max(lsp, axis=1) * 0.5
sos = [np.where(lsp[i,:] > half_max[i])[0][0] for i in range(lsp.shape[0])]
assert len(sos)==len(half_max)==lsp.shape[0]==len(npn)
npn.loc[:, 'lsp_sos'] = sos

# subset to unique NPN sites, getting site-means of mean_first_yes_doy
# NOTE: ideally we would have floristic data and weight each species by its
#       relative abundance here!
site_means = npn.groupby('site_id')['mean_first_yes_doy'].mean()
npn.loc[:, 'site_mean_first_yes_doy'] = [site_means[i] for i in npn['site_id']]

# plot relationship of our data to NPN and SIx data and compare to the
# relationship between them
fig = plt.figure(figsize=(12,4))
comp_cols = [['lsp_sos', 'site_mean_first_yes_doy',],
             ['lsp_sos', 'six_leaf_val'],
             ['site_mean_first_yes_doy', 'six_leaf_val'],
            ]
label_dict = {'lsp_sos': '$SOS_{NIR_{V}}$',
              'site_mean_first_yes_doy': 'NPN',
              'six_leaf_val': 'SI-x',
             }
for i in range(len(comp_cols)):
    ax = fig.add_subplot(1,3,i+1)
    # drop rows with missing SI-x values
    cols = comp_cols[i]
    keep_rows = np.all(pd.notnull(npn.loc[:, cols]), axis=1)
    subnpn = npn.loc[keep_rows, :]
    # build a regression; if comparison involves SIx data then include
    # a quadratic term, to account for the decreasing-increase form of the relxn
    X = subnpn.loc[:, cols[0]].values.reshape((len(subnpn), 1))
    quadratic = 'six_leaf_val' in cols
    if quadratic:
        X = np.hstack((X, (subnpn.loc[:, cols[0]].values**2).reshape(
                                                            (len(subnpn), 1))))
    y = subnpn.loc[:, cols[1]].values
    mod = sm.OLS(endog=y, exog=X, hasconst=False).fit()
    pred_X = np.vstack([np.arange(0, 365, 0.5)**i for i in range(1,3)]).T
    if not quadratic:
        pred_X = pred_X[:, 0].reshape((pred_X.shape[0], 1))
    pred_y = mod.predict(exog=pred_X)
    coeffs = mod.params
    coffs_pvals = mod.pvalues
    r2 = mod.rsquared
    pval = mod.f_pvalue
    ax.scatter(subnpn.loc[:, cols[0]],
               subnpn.loc[:, cols[1]],
               c='black',
               s=2,
               alpha=0.8,
              )
    ax.plot(pred_X[:, 0],
            pred_y,
            '-r',
            alpha=0.6,
            linewidth=0.75,
           )
    predictor = f"{np.round(coeffs[0], 2)} {label_dict[cols[0]]}"
    if quadratic:
        joiner = {True: '', False: '+'}[coeffs[1]<0]
        predictor = joiner.join([predictor, (f"{np.round(coeffs[1], 4)} "
                                            f"{label_dict[cols[0]]}")])
        predictor += "$^{2}$"
    eqxn_txt = f"{label_dict[cols[1]]} ~ {predictor}"
    r2_pval_txt = "$R^2 = $" f"{np.round(r2, 2)} (P={np.round(pval, 2)})"
    ax.text(10, 345, eqxn_txt, size=9)
    ax.text(10, 320, r2_pval_txt, color='red', size=9)
    ax.set_xlim(0, 365)
    ax.set_ylim(0, 365)
    ax.set_xlabel(label_dict[cols[0]], fontdict={'fontsize': 10})
    ax.set_ylabel(label_dict[cols[1]], fontdict={'fontsize': 10})
fig.subplots_adjust(left=0.06,
                    right=0.98,
                    bottom=0.12,
                    top=0.96,
                    wspace=0.25,
                   )


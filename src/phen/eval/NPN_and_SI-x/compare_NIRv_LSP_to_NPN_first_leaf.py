import os
import re
import sys
import pprint
import numpy as np
import pandas as pd
import rasterio as rio
import geopandas as gpd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from shapely.geometry import Point

# local modules
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
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

# subset to points where we have data
not_missing = np.all(pd.notnull(lsp), axis=1)
npn = npn.loc[not_missing, ]
lsp = lsp[not_missing, :]
assert lsp.shape[0] == len(npn)

print(("\n\nAfter accounting for masked LSP pixels, "
       f"data from {len(npn['site_id'].unique())} NPN "
       "sites remain.\n\n"))


# estimate start of season (SOS) as doy when LSP exceeds 50% of its amplitude
# (a la Melaas et al. 2016, but I take it this is standard practice)
half_max = np.percentile(lsp, 50, axis=1)
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
             ['six_leaf_val', 'site_mean_first_yes_doy'],
             ['lsp_sos', 'six_leaf_val'],
            ]
label_dict = {'lsp_sos': 'SOS$_{NIR_{V}}$',
              'site_mean_first_yes_doy': 'first leaf$_{NPN}$',
              'six_leaf_val': 'SOS$_{SI-x}$',
             }
letter_labs = ['A.', 'B.', 'C.']
for i in range(len(comp_cols)):
    ax = fig.add_subplot(1,3,i+1)
    # drop rows with missing SI-x values
    cols = comp_cols[i]
    print('*'*80)
    print('*'*80)
    print('*'*80)
    print(f"\n\nAnalyzing {cols[1]} ~ {cols[0]}...\n")
    keep_rows = np.all(pd.notnull(npn.loc[:, cols]), axis=1)
    subnpn = npn.loc[keep_rows, :]
    # simple linear regression
    X = subnpn.loc[:, cols[0]].values.reshape((-1, 1))
    X = sm.add_constant(X)
    y = subnpn.loc[:, cols[1]].values
    mod = sm.OLS(endog=y, exog=X).fit()
    print(mod.summary2())
    pred_X = np.arange(0, 365, 0.5).reshape((-1, 1))
    pred_X = np.hstack((np.ones(pred_X.shape), pred_X))
    pred_y = mod.predict(exog=pred_X)
    coeffs = mod.params
    r2 = mod.rsquared
    pval = mod.f_pvalue
    ax.scatter(subnpn.loc[:, cols[0]],
               subnpn.loc[:, cols[1]],
               c='black',
               s=1,
               alpha=0.8,
              )
    ax.plot(pred_X[:, 1],
            pred_y,
            '-r',
            alpha=0.6,
            linewidth=0.75,
           )
    predictor = f"{np.round(coeffs[1], 2)} {label_dict[cols[0]]}"
    intercept= f"{np.round(coeffs[0], 2)}"
    eqxn_txt = f"{label_dict[cols[1]]} ~\n{intercept} + {predictor}"
    r2_pval_txt = "$R^2 = $" + "%0.2f" % (np.round(r2, 2)) + f" (P≤{np.round(pval, 3)})"
    ax.text(175, 320, eqxn_txt, size=11)
    ax.text(175, 290, r2_pval_txt, color='red', size=11)
    ax.text(15, 330, letter_labs[i], size=24, weight='bold')
    ax.set_xlim(0, 365)
    ax.set_ylim(0, 365)
    ax.set_xlabel(label_dict[cols[0]], fontdict={'fontsize': 13})
    ax.set_ylabel(label_dict[cols[1]], fontdict={'fontsize': 13})
    ax.tick_params(labelsize=11)
fig.subplots_adjust(left=0.06,
                    right=0.98,
                    bottom=0.13,
                    top=0.96,
                    wspace=0.25,
                   )
fig.savefig(os.path.join(phf.FIGS_DIR,
                         'FIG_SUPP_NIRv_NPN_SI-x_SOS_comparison.png'), dpi=600)


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
from pykrige.uk import UniversalKriging

# local modules
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf
from MMRR import MMRR

#------------------------------------------------------------
# TODO:

# get rid of quadratic! (logistic instead? unnecessary?)

# what to do about very late 'first's?
'''
fig, ax = plt.subplots(1,1)
nam.plot(color='none', ax=ax)
npn[npn['mean_first_yes_doy']>300].plot('mean_first_yes_doy', ax=ax, zorder=1)
print(npn[npn['mean_first_yes_doy']>300].loc[:, ['genus', 'species']])
'''

# how/if to deal with spatial autocorrelation?

#------------------------------------------------------------


# data dir for saving large outputs
npn_data_dir = "/media/deth/SLAB/diss/3-phn/npn/"

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
        "of lsp consists of missing values"))

# subset to points where we have data
not_missing = np.all(pd.notnull(lsp), axis=1)
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
    # simple linear regression
    X = subnpn.loc[:, cols[0]].values.reshape((-1, 1))
    X = sm.add_constant(X)
    y = subnpn.loc[:, cols[1]].values
    mod = sm.OLS(endog=y, exog=X).fit()
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
    #edges = np.linspace(0, 365, int(365/5)+1)
    #heatmap, xedges, yedges = np.histogram2d(subnpn.loc[:, cols[0]],
    #                                         subnpn.loc[:, cols[1]],
    #                                         bins=edges,
    #                                        )
    #extent = [0, 365, 0, 365]
    #ax.imshow(heatmap.T, extent=extent, origin='lower', cmap='Greys')
    ax.plot(pred_X[:, 1],
            pred_y,
            '-r',
            alpha=0.6,
            linewidth=0.75,
           )
    predictor = f"{np.round(coeffs[1], 2)} {label_dict[cols[0]]}"
    intercept= f"{np.round(coeffs[0], 2)}"
    eqxn_txt = f"{label_dict[cols[1]]} ~ {intercept} + {predictor}"
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
fig.savefig('NIRv_NPN_SI-x_SOS_comparison.png', dpi=600)


# MMRR analysis

def calc_doy_diff(doy0, doy2):
    '''
    calculate the distance, in number of days, between 2 numericaly days of year
    '''
    # get the lesser of the distance back to the earlier day of year or
    # forward to the same day of year next year
    d = sorted([doy1, doy2])
    dist = np.min((d[1]-d[0], d[0] + 365 - d[1]))
    return dist

assert np.sum(pd.isnull(npn.loc[:, ['lsp_sos',
                                    'site_mean_first_yes_doy']].values)) == 0

# set the numpy.random seed
np.random.seed(1)

# read dist matrices, if saved:
lsp_fn = os.path.join(npn_data_dir, f"dist_mat_LSP_SOS.txt")
npn_fn = os.path.join(npn_data_dir, f"dist_mat_NPN_leaves.txt")
if os.path.isfile(lsp_fn) and os.path.isfile(npn_fn):
    lsp_sos_dists = np.loadtxt(lsp_fn)
    npn_flf_dists = np.loadtxt(npn_fn)
else:
    # calculate first-leaf date distances for both LSP and NPN datasets
    # (i.e., numbers of calendar days separating the doy values)
    print('\n\tcalculating first-leaf distances...')
    lsp_sos_dists = np.ones([len(npn)]*2)*np.nan
    npn_flf_dists = np.ones([len(npn)]*2)*np.nan
    for i in range(len(npn)):
        lsp_sos_i = npn.iloc[i, :]['lsp_sos']
        npn_flf_i = npn.iloc[i, :]['site_mean_first_yes_doy']
        for j in range(len(npn)):
            if i==j:
                lsp_sos_dists[i, j] = 0
                npn_flf_dists[i, j] = 0
            else:
                lsp_sos_j = npn.iloc[j, :]['lsp_sos']
                lsp_sos_dist = calc_doy_diff(lsp_sos_i, lsp_sos_j)
                lsp_sos_dists[i, j] = lsp_sos_dists[j, i] = lsp_sos_dist
                npn_flf_j = npn.iloc[j, :]['site_mean_first_yes_doy']
                npn_flf_dist = calc_doy_diff(npn_flf_i, npn_flf_j)
                npn_flf_dists[i, j] = npn_flf_dists[j, i] = npn_flf_dist
    assert np.sum(pd.isnull(lsp_sos_dists)) == 0
    assert np.sum(pd.isnull(npn_flf_dists)) == 0
    assert np.all(lsp_sos_dists == lsp_sos_dists.T)
    assert np.all(npn_flf_dists == npn_flf_dists.T)
    for label, dist_mat in zip(['LSP_SOS', 'NPN_leaves'],
                               [lsp_sos_dists, npn_flf_dists],):
        np.savetxt(os.path.join(npn_data_dir,
                                f"dist_mat_{label}.txt"),
                   dist_mat)

# run MMRR and store results
# run the MMRR model and print results
print('\n\trunning MMRR model...')
res = MMRR(Y=npn_flf_dists,
           X=[lsp_sos_dists],
           Xnames=['lsp_sos_dist'],
           # NOTE: MMRR will standardize lower-triangular distance values, and thus
           #       returns coefficient values as beta-coefficients
           standardize=True,
           intercept=True,
           nperm=999,
          )

# save MMRR results
pd.DataFrame(res, index=[0]).to_csv('NPN_MMRR_results.csv', index=False)

# plot MMRR results
lsp_sos_vals = lsp_sos_dists[np.tril_indices(lsp_sos_dists.shape[0], k=-1)]
npn_flf_vals = npn_flf_dists[np.tril_indices(npn_flf_dists.shape[0], k=-1)]
for a in [lsp_sos_vals, npn_flf_vals]:
    assert a.size == (len(npn) * (len(npn)-1))/2
    assert np.sum(pd.isnull(a)) == 0
    assert np.min(a) == 0
    assert np.max(a) <= 365/2

fig, ax = plt.subplots(1,1, figsize=(5,5))
# simple linear regression
X = lsp_sos_vals.ravel().reshape((-1, 1))
X = sm.add_constant(X)
y = npn_flf_vals.ravel().reshape((-1, 1))
mod = sm.OLS(endog=y, exog=X).fit()
pred_X = np.arange(0, 365/2, 0.1).reshape((-1, 1))
pred_X = sm.add_constant(pred_X)
pred_y = mod.predict(exog=pred_X)
coeffs = mod.params
r2 = mod.rsquared
pval = mod.f_pvalue
#ax.scatter(lsp_sos_vals.ravel(),
#           npn_flf_vals.ravel(),
#           c='black',
#           s=1,
#           alpha=0.8,
#          )
edges = np.linspace(0, 365/2, int(365/2/5)+1)
heatmap, xedges, yedges = np.histogram2d(lsp_sos_vals.ravel(),
                                         npn_flf_vals.ravel(),
                                         bins=edges,
                                        )
extent = [0, 365/2, 0, 365/2]
ax.imshow(heatmap.T, extent=extent, origin='lower', cmap='Greys')
ax.plot(pred_X[:, 1],
        pred_y,
        '-r',
        alpha=0.6,
        linewidth=0.75,
       )
predictor = f"{np.round(coeffs[1], 2)} "
predictor = predictor + "$SOS_{NIR_{V}}"
intercept= f"{np.round(coeffs[0], 2)}"
eqxn_txt = f"{label_dict[cols[1]]} ~ {intercept} + {predictor}"
r2_pval_txt = "$R^2 = $" f"{np.round(r2, 2)} (P={np.round(pval, 2)})"
ax.text(10, 50, eqxn_txt, size=9)
ax.text(10, 45, r2_pval_txt, color='red', size=9)
ax.set_xlim(0, 365/2)
ax.set_ylim(0, 365/2)
ax.set_xlabel('$SOS_{NIR_{V}}$ distance (days)', fontdict={'fontsize': 10})
ax.set_ylabel('NPN first-leaf distance (days)', fontdict={'fontsize': 10})
fig.subplots_adjust(left=0.06,
                    right=0.98,
                    bottom=0.12,
                    top=0.96,
                    wspace=0.25,
                   )
fig.savefig('NPN_MMRR_results.png', dpi=600)


# krig and display results
npn_krig = npn[npn['longitude'].between(-130, -60)]
# TODO: DELETE NEXT ROW TO RUN FOR FULL DATA
npn_krig = npn_krig.iloc[::10, :]
krig = UniversalKriging(npn_krig.to_crs(8857)['longitude'],
                        npn_krig.to_crs(8857)['latitude'],
                        npn_krig['site_mean_first_yes_doy'],
                       )
minx, miny = np.min(npn_krig.to_crs(8857).bounds.loc[:, ['minx', 'miny']],
                    axis=0).values
maxx, maxy = np.max(npn_krig.to_crs(8857).bounds.loc[:, ['maxx', 'maxy']],
                    axis=0).values
# krig values to a ~5km grid
gridx = np.arange(minx, maxx, 5000)
gridy = np.arange(miny, maxy, 5000)
krig_grid, krig_ss = krig.execute("grid", gridx, gridy)
fig, ax = plt.subplots(1,1, figsize=(5,5))
ax.imshow(krig_grid, cmap='twilight', zorder=0)
world.to_crs(8857).plot(color='none', edgecolor='black', alpha=0.5, zorder=1)

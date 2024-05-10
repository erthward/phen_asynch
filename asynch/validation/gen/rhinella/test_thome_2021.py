"""
test our fitted LSP patterns against the data from
Thomé et al. 2021 (DOI 10.1038/s41437-021-00460-7)
by redoing their MMRR test but using our LSP
patterns to test the Asynchrony of Seasons hypothesis
in lieu of their monthly climatic LDA variables
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
import rasterio as rio
import rioxarray as rxr
import matplotlib.pyplot as plt
import sys
from sklearn.cluster import KMeans


sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf
from MMRR import MMRR


# set the numpy.random seed
np.random.seed(1)


# whether to interpolate missing LSP coefficients, 
# and if so, how far away can rasterio look for neighbor to interpolate from?
# coefficient raster?
interp_sea_data = False
neigh_dist_sea_fill_tol = 5

# how many MMRR permutations to use
MMRR_nperm = 999


# read genetic distance matrix
print('\nreading genetic data...')
# NOTE: simple Euclidean genetic distance matrix calculated using adegenet in R
gen_dist = pd.read_csv('./dryad_archive/2-Mantel_test/Rgranulosa_Mantel.csv')
# reformat row and col names to reflect just the voucher code
# (without the lingering haplotype integer)
gen_dist.index = [i.split('_')[0] for i in gen_dist.index]
gen_dist.columns = [c.split('_')[0] for c in gen_dist.columns]
# check voucher IDs on rows and columns are unique and identical
assert len(np.unique(gen_dist.index)) == len(gen_dist)
assert len(np.unique(gen_dist.columns)) == len(gen_dist.columns)
assert np.all(np.array([*gen_dist.index]) == np.array([*gen_dist.columns]))
# check diagonal is all zeros
assert np.all([gen_dist.iloc[i, i] == 0 for i in range(gen_dist.shape[0])])

# read CSV with sample voucher codes and locations
print('\nreading geographic data...')
geo = pd.read_csv('./41437_2021_460_MOESM1_ESM_TABLES1.csv')
# get rid of white space in column names and values
geo.columns = [c.strip() for c in geo.columns]
# rename lon and lat
geo.rename(mapper={'Longitude': 'lon', 'Latitude': 'lat'}, axis=1, inplace=True)
geo['voucher'] = [v.strip() for v in geo['voucher'].values]
# NOTE: MTCT1468 in gen_dist rows/columns must match
#       MTCT1368 in geo['voucher'], as they're the only
#       ones that reciprocally have no match in each other.
#       FIX THAT!
geo.loc[geo['voucher'] == 'MTCT1368', 'voucher'] = 'MTCT1468'
# check we have matching voucher IDs in both datasets
assert set(geo['voucher'].values) == set(gen_dist.index)


# set up the geographic, environmental, and seasonal distance matrices
geo_dist = np.ones([geo.shape[0]]*2) * np.nan

# loop over voucher codes in gen_dist, to ensure that rows & cols
# match the rows & cols in gen_dist
print('\ncalculating geographic distances...')
g = pyproj.Geod(ellps='WGS84')
pts = []
for i, vouch_i in enumerate(gen_dist.index):
    lon_i, lat_i = geo.loc[geo['voucher'] == vouch_i, ['lon', 'lat']].values[0]
    pts.append([lon_i, lat_i])
    for j, vouch_j in enumerate(gen_dist.columns):
        # calculate geographic distance
        if i == j:
            dist_ij = 0
        else:
            lon_j, lat_j = geo.loc[geo['voucher'] == vouch_j, ['lon', 'lat']].values[0]
            (az_ij, az_ji, dist_ij) = g.inv(lon_i, lat_i, lon_j, lat_j)
        geo_dist[i, j] = dist_ij
        geo_dist[j, i] = dist_ij


pts = np.array(pts)
# calculate environmental distance for just the 4 variables they list
# (BIO1, BIO4, BIO12, BIO15)
print('\ncalculating environmental distances...')
nodata_val = rio.open(phf.BIOCLIM_INFILEPATHS[0]).nodata
env_dist = phf.calc_pw_clim_dist_mat(pts,
                                     nodata_val=nodata_val,
                                     vars_list=['bio_1',
                                                'bio_4',
                                                'bio_12',
                                                'bio_15',
                                               ],
                                    )

# calculate seasonal distance
# NOTE: needs to be based on standardized time series,
#       to focus on timing rather than fitted NIRV values
print('\ncalculating seasonal distances...')
sea_dist = phf.get_raster_info_points(phf.COEFFS_STRICT_FILE,
                                      pts,
                                      'ts_pdm',
                                      standardize=True,
                                      fill_nans=interp_sea_data,
                                      fill_tol=neigh_dist_sea_fill_tol,
                                     )
print((f"\n\n{np.round(100*np.mean(pd.isnull(sea_dist)), 2)} % "
        "of sea_dist consists of missing values"))

# check all diagonals are zeros
assert np.all([geo_dist[i, i] == 0 for i in range(geo_dist.shape[0])])
assert np.all([env_dist[i, i] == 0 for i in range(env_dist.shape[0])])
assert np.all([sea_dist[i, i] == 0 for i in range(sea_dist.shape[0])])

# drop rows and cols in gen_dist, geo_dist, and env_dist
# if they feature NaNs in sea_dist
# (happens because our raster has NaNs where pixels were filtered out)
print('\ndropping rows with missing LSP data...')
missing_sea = (np.sum(pd.isnull(sea_dist), axis=0) == sea_dist.shape[0] - 1)
missing_sea_vouch = gen_dist.index[np.where(missing_sea)].values
print(("\n\nThe following voucher IDs occur at locations where seasonality "
       "data is missing and was not interpolated:\t"
       f"{'   '.join(missing_sea_vouch)}\n\n"))
fig, ax = plt.subplots(1,1)
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
world.plot(color='none', ax=ax)
missing_sea_color = ['blue' if v not in missing_sea_vouch else 'red' for
                                                    v in geo['voucher'].values]
geo['missing_sea'] = missing_sea_color
ax.scatter(x=geo['lon'], y=geo['lat'], c=geo['missing_sea'], s=25, alpha=0.5)
ax.set_title('red points missing seasonality data in our LSP dataset')
ax.set_xlim(-47.5, -32.5)
ax.set_ylim(-25, 0)
fig.show()

gen_dist = gen_dist.values
gen_dist = gen_dist[~missing_sea].T[~missing_sea]
geo_dist = geo_dist[~missing_sea].T[~missing_sea]
env_dist = env_dist[~missing_sea].T[~missing_sea]
sea_dist = sea_dist[~missing_sea].T[~missing_sea]
for i, mat in enumerate([gen_dist, geo_dist, env_dist, sea_dist]):
    assert np.all(mat == mat.T), f"matrix {i} failed!"
    assert np.all(mat.shape == gen_dist.shape), f"matrix {i} failed!"


# TODO: need to include resistance/geo barrier, given that they provide no
#       data/code for it and the var was insignificant in their model anyhow?



# TODO: need to include instability, given that they provide no
#       data/code for it and the var was insignificant in their model anyhow?



# run the MMRR model and print results
print('\nrunning MMRR model...')
res = MMRR(Y=gen_dist,
           X=[geo_dist, env_dist, sea_dist],
           Xnames=['geo_dist', 'env_dist', 'sea_dist'],
           # NOTE: MMRR will standardize lower-triangular distance values, and thus
           #       returns coefficient values as beta-coefficients
           standardize=True,
           intercept=True,
           nperm=MMRR_nperm,
          )

print(res)


# map genetic sampling locations on left
# and plot fitted seasonal time series on right,
# in both instances colored by a K-means clustering
# with K=2, using colors taken from Thomé et al. plots
sea_ts = phf.get_raster_info_points(phf.COEFFS_STRICT_FILE,
                                    pts[~missing_sea,:],
                                    'ts',
                                    standardize=True,
                                    fill_nans=interp_sea_data,
                                    fill_tol=neigh_dist_sea_fill_tol,
                                   )
km = KMeans(n_clusters=2).fit(sea_ts)
fig, axs = plt.subplots(1,2)
red = '#ca1957'
blue = '#2d5098'
colors = np.array([blue, red])
# plot map
ax = axs[0]
world.plot(color='none', zorder=1, ax=ax)
ax.scatter(pts[~missing_sea, 0],
           pts[~missing_sea, 1],
           c=colors[km.labels_],
           zorder=2,
          )
ax.set_xlim(-47.5, -32.5)
ax.set_ylim(-25, 0)
ax.set_xticks(())
ax.set_yticks(())
ax.set_title('sampling points,\nclustered by LSP pattern',
             fontdict={'fontsize': 19})
# plot LSP time series
ax = axs[1]
for i in range(sea_ts.shape[0]):
    ax.plot(sea_ts[i, :], color=colors[km.labels_[i]], alpha=0.7)
xax_ticks = [0, 90, 180, 271, 364]
xax_ticklabs = ['Jan', 'Apr', 'Jul', 'Oct', 'Jan']
ax.set_xticks(xax_ticks)
ax.set_xticklabels(xax_ticklabs, size=14)
ax.set_yticklabels(ax.get_yticklabels(), size=14)
ax.set_xlim(0, 364)
ax.set_xlabel('day of calendar year', fontdict={'fontsize': 19})
ax.set_ylabel('normalized LSP value', fontdict={'fontsize': 19})
ax.set_title('LSP patterns at sampling points,\ncolored by cluster',
             fontdict={'fontsize': 19})
fig.show()


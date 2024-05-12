"""
test our fitted LSP patterns against all
GenBank Xiphorhynchus fuscus data used
in the only eastern Brazilian analysis in
Qunitero et al. 2014 (DOI 10.1086/677261),
originally from Cabanne et al. 2008
(DOI 10.1016/j.ympev.2008.09.013)
by running an MMRR and using our LSP
patterns to test the Asynchrony of Seasons Hypothesis
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
import rasterio as rio
import rioxarray as rxr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import re
import pprint
from sklearn.cluster import KMeans
from Bio import SeqIO


sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf
from MMRR import MMRR


print('\n\nRunning LSP MMRR test for Xiphorhynchus fuscus...')


# set the numpy.random seed
np.random.seed(1)

# whether to interpolate missing LSP coefficients, 
# and if so, how far away can rasterio look for neighbor to interpolate from?
# coefficient raster?
interp_sea_data = False
neigh_dist_sea_fill_tol = 5

# how many MMRR permutations to use
MMRR_nperm = 999


# load country boundaries
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                       'NewWorldFile_2020.shp'))

# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR,f)))
subnational = pd.concat(subnational)


# read genetic distance matrix
print('\nreading genetic data...')
gen_dist = pd.read_csv('./all_xiphorhynchus_fuscus_dist.csv')
# reformat row and col names to reflect just the voucher code
gen_dist.index = [i.split('.')[0] for i in gen_dist.index]
gen_dist.columns = [c.split('.')[0] for c in gen_dist.columns]
# check voucher IDs on rows and columns are unique and identical
assert len(np.unique(gen_dist.index)) == len(gen_dist)
assert len(np.unique(gen_dist.columns)) == len(gen_dist.columns)
assert np.all(np.array([*gen_dist.index]) == np.array([*gen_dist.columns]))
# check diagonal is all zeros
assert np.all([gen_dist.iloc[i, i] == 0 for i in range(gen_dist.shape[0])])

# read CSV with sample voucher codes and locations
print('\nreading geographic data...')
geo = pd.read_csv('./all_xiphorhynchus_fuscus_locs.csv')
# calculate decimal degree columns
# NOTE: only missing vals are geo coord secs in most rows, so replace w/ 0s
geo = geo.replace(np.nan, 0)
# NOTE: multiply all by -1 because all in S and W hemispheres
geo['lat'] = -1*(geo['lat_d'] + ((geo['lat_m'] + (geo['lat_s']/60))/60))
geo['lon'] = -1*(geo['lon_d'] + ((geo['lon_m'] + (geo['lon_s']/60))/60))
# check we have matching voucher IDs in both datasets
assert set(geo['id'].values) == set(gen_dist.index)

# set up the geographic, environmental, and seasonal distance matrices
geo_dist = np.ones([geo.shape[0]]*2) * np.nan

# loop over voucher codes in gen_dist, to ensure that rows & cols
# match the rows & cols in gen_dist
print('\ncalculating geographic distances...')
g = pyproj.Geod(ellps='WGS84')
pts = []
for i, vouch_i in enumerate(gen_dist.index):
    lon_i, lat_i = geo.loc[geo['id'] == vouch_i, ['lon', 'lat']].values[0]
    pts.append([lon_i, lat_i])
    for j, vouch_j in enumerate(gen_dist.columns):
        # calculate geographic distance
        if i == j:
            dist_ij = 0
        else:
            lon_j, lat_j = geo.loc[geo['id'] == vouch_j, ['lon', 'lat']].values[0]
            (az_ij, az_ji, dist_ij) = g.inv(lon_i, lat_i, lon_j, lat_j)
        geo_dist[i, j] = dist_ij
        geo_dist[j, i] = dist_ij


pts = np.array(pts)
# calculate environmental distance across all Worldclim vals
print('\ncalculating environmental distances...')
nodata_val = rio.open(phf.BIOCLIM_INFILEPATHS[0]).nodata
env_dist = phf.calc_pw_clim_dist_mat(pts,
                                     nodata_val=nodata_val,
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
countries.plot(color='none', ax=ax)
missing_sea_color = ['blue' if v not in missing_sea_vouch else 'red' for
                                                    v in geo['id'].values]
geo['missing_sea'] = missing_sea_color
ax.scatter(x=geo['lon'], y=geo['lat'], c=geo['missing_sea'], s=25, alpha=0.5)
ax.set_title('red points missing seasonality data in our LSP dataset')
ax.set_xlim(-57.5, -34)
ax.set_ylim(-30, 0)
fig.savefig('xiphorhynchus_fuscus_sites_missing_LSP_data.png', dpi=300)

gen_dist = gen_dist.values
gen_dist = gen_dist[~missing_sea].T[~missing_sea]
geo_dist = geo_dist[~missing_sea].T[~missing_sea]
env_dist = env_dist[~missing_sea].T[~missing_sea]
sea_dist = sea_dist[~missing_sea].T[~missing_sea]
for i, mat in enumerate([gen_dist, geo_dist, env_dist, sea_dist]):
    assert np.all(mat == mat.T), f"matrix {i} failed!"
    assert np.all(mat.shape == gen_dist.shape), f"matrix {i} failed!"


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
print('\nMMRR results:\n')
pprint.pprint(res)


# map genetic sampling locations on left
# and plot fitted seasonal time series on right,
# all colored by a K-means clustering
# with K=2, using colors taken from Thomé et al. plots
sea_ts = phf.get_raster_info_points(phf.COEFFS_STRICT_FILE,
                                    pts[~missing_sea,:],
                                    'ts',
                                    standardize=True,
                                    fill_nans=interp_sea_data,
                                    fill_tol=neigh_dist_sea_fill_tol,
                                   )
km_lsp = KMeans(n_clusters=2).fit(sea_ts)
# flip the two clusters' labels so that the first point is cluster 0
# (labels are arbitrary anyhow, but this aligns the colors with Thomé Fig. 5,
# i.e. blue in North)
if km_lsp.labels_[0] == 1:
    km_lsp.labels_ = np.int16(km_lsp.labels_ == 0)
fig  = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(2,2)
red = '#ca1957'
blue = '#2d5098'
colors = np.array([blue, red])
# plot map
ax = plt.subplot(gs[:, 0])
subnational.plot(color='none',
                 edgecolor='black',
                 zorder=0,
                 ax=ax,
                 alpha=0.6,
                 linewidth=0.05,
                )
countries.plot(color='none',
                            edgecolor='black',
                            linewidth=0.2,
                            zorder=1,
                            ax=ax,
                            )
ax.scatter(pts[~missing_sea, 0],
           pts[~missing_sea, 1],
           s=11,
           c=colors[km_lsp.labels_],
           edgecolor='black',
           linewidth=0.5,
           alpha=0.6,
           zorder=2,
          )
ax.set_xlim(-57.5, -34)
ax.set_ylim(-30, 0)
ax.set_xticks(())
ax.set_yticks(())
ax.set_title('sampling points,\nclustered by LSP pattern',
             fontdict={'fontsize': 10})
# plot LSP time series clusters
ax = plt.subplot(gs[0, 1])
for i in range(sea_ts.shape[0]):
    ax.plot(sea_ts[i, :],
            color=colors[km_lsp.labels_[i]],
            alpha=0.7,
            linewidth=0.5,
           )
xax_ticks = [0, 90, 180, 271, 364]
xax_ticklabs = ['Jan', 'Apr', 'Jul', 'Oct', 'Jan']
ax.set_xticks(xax_ticks)
ax.set_xticklabels(xax_ticklabs, size=7)
assert np.round(np.max(sea_ts), 0) == 2
assert np.round(np.min(sea_ts), 0) == -2
ax.set_yticks([*range(-2, 3)],
              [str(n) for n in range(-2, 3)],
              size=7,
             )
ax.set_xlim(0, 364)
ax.set_xlabel('day of calendar year', fontdict={'fontsize': 8})
ax.set_ylabel('normalized LSP value', fontdict={'fontsize': 8})
ax.set_title('colored by LSP cluster',
             fontdict={'fontsize': 10})

# load genetic data
seqs = []
for rec in SeqIO.parse('./all_xiphorhynchus_fuscus_ALIGNED.fasta',
                       format='fasta'):
    seqs.append(str(rec.seq))
assert len({*[len(s) for s in seqs]}) == 1
# delete sites with any missing bases
miss = []
for i in range(len(seqs[0])):
    for s in seqs:
        if s[i] in ['-', 'n']:
            miss.append(i)
miss = list({*miss})
keep = [i for i in range(len(seqs[0])) if i not in miss]
seqs = np.array([np.array([*s])[keep] for s in seqs])
# drop vouchers with missing LSP data
seqs = seqs[~missing_sea, :]
# remap bases to numbers, for clustering
seqs_n = np.ones(seqs.shape)*np.nan
assert np.all(sorted(np.unique(seqs)) == np.array(['a', 'c', 'g', 't']))
for n, base in enumerate(['a', 'c', 'g', 't']):
    seqs_n[seqs==base] = n
# plot LSP time series, colored by genetic clusters
ax = plt.subplot(gs[1, 1])
km_gen = KMeans(n_clusters=2).fit(seqs_n)
# flip the two clusters' labels so that the first point is cluster 0
# (labels are arbitrary anyhow, but this aligns the colors with Thomé Fig. 5,
# i.e. blue in North)
if km_gen.labels_[0] == 1:
    km_gen.labels_ = np.int16(km_gen.labels_ == 0)
for i in range(sea_ts.shape[0]):
    ax.plot(sea_ts[i, :],
            color=colors[km_gen.labels_[i]],
            alpha=0.7,
            linewidth=0.5,
           )
xax_ticks = [0, 90, 180, 271, 364]
xax_ticklabs = ['Jan', 'Apr', 'Jul', 'Oct', 'Jan']
ax.set_xticks(xax_ticks)
ax.set_xticklabels(xax_ticklabs, size=7)
assert np.round(np.max(sea_ts), 0) == 2
assert np.round(np.min(sea_ts), 0) == -2
ax.set_yticks([*range(-2, 3)],
              [str(n) for n in range(-2, 3)],
              size=7,
             )
ax.set_xlim(0, 364)
ax.set_xlabel('day of calendar year', fontdict={'fontsize': 8})
ax.set_ylabel('normalized LSP value', fontdict={'fontsize': 8})
ax.set_title('colored by genetic cluster',
             fontdict={'fontsize': 10})
fig.subplots_adjust(hspace=0.5,
                    wspace=0.1,
                    left=0.0,
                    right=0.97,
                    bottom=0.1,
                    top=0.9,
                   )
fig.savefig('xiphorhynchus_fuscus_results.png', dpi=600)


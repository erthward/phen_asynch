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
import matplotlib.gridspec as gridspec
import os
import re
import sys
import pprint
from sklearn.cluster import KMeans


# TODO:
    # HOW TO DEAL WITH FACT THAT NOT A SINGLE LOCUS HAS DATA FOR ALL VOUCHERS??



sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf
from MMRR import MMRR


print('\n\nRunning LSP MMRR test for Rhinella granulosa...')


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

# loop over voucher codes, ordered by gen_dist rows,
# to ensure that rows & cols match the rows & cols in gen_dist
vouchers = np.array([*gen_dist.index])
print('\ncalculating geographic distances...')
g = pyproj.Geod(ellps='WGS84')
pts = []
for i, vouch_i in enumerate(vouchers):
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
missing_sea_vouch = vouchers[np.where(missing_sea)]
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
ax.set_xlim(-57.5, -34)
ax.set_ylim(-30, 0)
fig.savefig('rhinella_granulosa_sites_missing_LSP_data.png', dpi=300)

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
fig  = plt.figure(figsize=(9,5))
gs = gridspec.GridSpec(2,2)
red = '#ca1957'
blue = '#2d5098'
colors = np.array([blue, red])
# plot map
ax = plt.subplot(gs[:, 0])
phf.plot_juris_bounds(ax,
                      lev1_zorder=0,
                      lev1_alpha=0.6,
                      lev1_linewidth=0.05,
                      lev0_zorder=1,
                      lev0_alpha=1,
                      lev0_linewidth=0.2,
                      strip_axes=False,
                      crs=4326,
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

## load genetic data
#with open('./dryad_archive/1-SNP_loci_pyrad/Rgranulosa.loci', 'r') as f:
#    all_loci_txt = f.read().split('|\n')
#    all_loci = {}
#    for loc_txt in all_loci_txt:
#        print(len(loc_txt.split('\n')))
#        if re.search("(?<=\|)\d+", loc_txt) is not None:
#            loc_num = int(re.search("(?<=\|)\d+", loc_txt).group())
#            loc_txt_split = loc_txt.split('\n//')
#            assert len(loc_txt_split) == 2
#            loc_seqs = {}
#            for seq in loc_txt_split[0].split('\n'):
#                seq_split = seq.split()
#                assert len(seq_split) == 2
#                vouch = re.search("(?<=>)[A-Z0-9]*", seq_split[0]).group()
#                seq_str = seq_split[1].strip()
#                loc_seqs[vouch] = seq_str
#            all_loci[loc_num] = loc_seqs
#
## process genetic data
#concat_seqs = {vouch: [] for vouch in vouchers}
#base_num_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
#for loc_num, seqs in all_loci.items():
#    assert len({*[len(seq) for seq in seqs.values()]}) == 1
#    # delete sites with any missing bases
#    miss = []
#    for vouch in vouchers:
#        try:
#            seq = seqs[vouch]
#            for i, s in enumeate(seq):
#                if s[i] not in ['A', 'C' 'G', 'T']:
#                    miss.append(i)
#        except Exception:
#            pass
#    miss = list({*miss})
#    keep = [i for i in range(len([*seqs.values()][0])) if i not in miss]
#    for vouch in vouchers:
#        try:
#            filtered_seq = np.array([*seqs[vouch]])[keep]
#            num_seq = np.array([base_num_dict[b] for b in filtered_seq])
#            concat_seqs[vouch].append(num_seq)
#        except Exception:
#            concat_seqs[vouch].append(np.array([np.nan]*len(keep)))
#concat_seqs = {vouch: np.concatenate(seqs) for vouch, seqs in concat_seqs.items()}
#assert len(np.unique([len(seq) for seq in concat_seqs.values()])) == 1
## keep only loci with data for all vouchers
#gen_dat = np.ones((len(vouchers), len([*concat_seqs.values()][0])))*np.nan
#for i, vouch in enumerate(vouchers):
#    gen_dat[i, :] = concat_seqs[vouch]
#keep_cols = pd.isnull(gen_dat).sum(axis=0) == 0
#gen_dat = gen_dat[:, keep_cols]
## TODO: WHAT TO DO NOW THAT NO LOCI REMAIN!?

# plot LSP time series, colored by genetic clusters
ax = plt.subplot(gs[1, 1])
#km_gen = KMeans(n_clusters=2).fit(gen_dat)
km_gen = KMeans(n_clusters=2).fit(gen_dist)
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
fig.savefig('rhinella_granulosa_results.png', dpi=600)


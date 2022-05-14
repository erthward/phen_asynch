#!/usr/bin/env python
# coding: utf-8


# TODO:

    # do I need to restrict to only clusters whose polygons cover a certain
    # area, or contain a certain minimum number of points?

    # decide whether or not I should be standardizing seasonal time series
    # (or just rerun to compare)


# py packages
import geopandas as gpd
import numpy as np
import glob
import json
import time
import os
#from osgeo import gdal
import random
import pyproj
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
from copy import deepcopy
from scipy.spatial import KDTree
import seaborn as sns
import math
import itertools
import alphashape
import rasterio as rio
import rioxarray as rxr
import matplotlib.pyplot as plt
from descartes import PolygonPatch
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn import cluster
from shapely.ops import unary_union
from shapely.geometry import Polygon, Point, MultiPolygon
from matplotlib.patches import Polygon as mplPolygon
from scipy import stats
from scipy.spatial import ConvexHull
from scipy.spatial import distance
from collections import Counter as C

# local modules
import helper_fns as hf
from MMRR import MMRR


# CLUSTER PARAMS:
# dbscan params
dbscan_eps = 2
dbscan_min_samples_frac_all = 0.1# alpha hull param
alpha = 1.25

# ANALYSIS PARAMS:
# a number <= number of random points desired
n_pts = 1000
# standardize the seasonal time series before calculating their distances?
standardize_ts = True
# file suffix depending on standardization
if standardize_ts:
    file_suffix = '_STANDARDIZED'
else:
    file_suffix = '_UNSTANDARDIZED'

# run MMRR tests?
run_MMRR = True
# make time series figure?
plot_ts = False
# make the boxplot?
make_boxplot = False


#####################################
# PART I: CLUSTER HIGH-ASYNCH REGIONS
#####################################

# load the country boundaries (to use for random point generation)
cntry = gpd.read_file(hf.COUNTRIES_DATA_DIR + 'countries.shp')

# load asynch data
asynch = rxr.open_rasterio('../../results/maps/NIRv_global_asynch.tif')[2]

# mask everything below Nth percentile
pctile = np.nanpercentile(asynch, 95)
maxes = asynch.where(asynch>=pctile, np.nan).clip(min=1, max=1)

# extract as points
X, Y = np.meshgrid(maxes.x.values, maxes.y.values)
coords = np.array([*zip(X.ravel(), Y.ravel())])
coords = coords[np.where(pd.notnull(maxes.values.ravel()))[0]]

# try out clustering algo
# NOTE: DBSCAN SEEMED LIKE THE BEST FIT OFF THE BAT,
#       BASED ON DEPICTION HERE:
#           https://scikit-learn.org/stable/modules/clustering.html#overview-of-clustering-methods
#       AND DESCRIPTION HERE:
#           https://scikit-learn.org/stable/modules/clustering.html#dbscan
db = cluster.DBSCAN(eps=dbscan_eps,
        min_samples=dbscan_min_samples_frac_all * maxes.shape[0]).fit(coords)

# NOTE: ADD 1 TO CLUSTER LABELS, SO THAT FIRST CLUSTER == 1, NOT 0
db.labels_ += 1

# get clusters' mean latitudes (to be able to plot in order, N to S
clust_mean_lats = {}
for clust_i in np.unique(db.labels_):
    if clust_i > 0:
        mean_lat = np.mean(coords[db.labels_==clust_i,1])
        clust_mean_lats[clust_i] = mean_lat

# remap cluster labels so that they go in latitudinal order (N->S)
NtoS_clust_labels = sorted(clust_mean_lats.keys(),
                           key=lambda k: clust_mean_lats[k])[::-1]
new_labels_dict = {0:0}
for i, new_label in enumerate(NtoS_clust_labels):
    new_labels_dict[new_label] = i+1
new_labels_list = []
for label in db.labels_:
    new_labels_list.append(new_labels_dict[label])

fig_map = plt.figure(figsize=(20,8))
gs = fig_map.add_gridspec(10, 2, width_ratios=[1,0.25])
if make_boxplot:
    ax_map = fig_map.add_subplot(gs[:,1])
else:
    ax_map = fig_map.add_subplot(gs[:,:])
cntry.to_crs(4326).plot(color='none', edgecolor='k', linewidth=0.25, ax=ax_map)
ax_map.scatter(coords[:, 0], coords[:, 1], c=db.labels_, cmap='tab20', s=3)
ax_map.scatter(coords[:, 0], coords[:, 1], c=db.labels_>0, cmap='Greys_r',
               alpha=db.labels_==0, s=3)

# get list of regions (i.e., cluster polygons)
regs = []
for clust_i in np.unique(db.labels_):
    # ignore noisy samples
    if clust_i > 0:
        print('\n\nGETTING POLYGON AND MAPPING FOR CLUSTER %i...\n\n' % clust_i)
        # indices of this cluster's samples
        inds = np.where(db.labels_==clust_i)[0]
        # inds = [*set([*np.where(db.labels_==clust_i)[0]]).intersection(
                 #set(db.core_sample_indices_))]
        # coords of this cluster's samples
        clust_coords = coords[inds,:]
        # get alpha-hull coords
        alpha_hull = alphashape.alphashape(clust_coords, alpha=alpha)
        regs.append(alpha_hull)
        # add to map
        ax_map.add_patch(PolygonPatch(alpha_hull, alpha=0.2, color='black'))
        cent = alpha_hull.centroid.coords.xy
        ax_map.text(cent[0][0], cent[1][0], str(clust_i), color='black',
                    size=14, weight='bold', alpha=0.85)
# coerce all to MultiPolygons (to be able to use for raster clipping
regs = [reg if isinstance(reg, MultiPolygon) else MultiPolygon([reg]) for reg in regs]


####################################
# PART II: DRAW PTS AND RUN ANALYSIS
####################################

# figure to plot ts for each reg separately
if plot_ts:
    fig_ts = plt.figure(figsize=(15,10))

# read in the Köppen data
kopp = rxr.open_rasterio('./clim_dist/Beck_rewrite_0p0083_5KM_AGG.tif')[0]
# mask out oceans
kopp = kopp.where(kopp>0, np.nan)

# draw random points for each region
regs_pts = [hf.generate_random_points_in_polygon(n_pts, reg) for reg in regs]

# create the figure
fig, ax = plt.subplots(1,1)
fig.suptitle(('seasonal distance vs. climatic distance, '
              'across latitude and asynchrony'), size=20)

# list of region names
reg_names = [str(i+1) for i in range(len(regs))]

# create colors for the 4 regions
reg_cols = [plt.cm.tab20(v/len(regs)) for v in range(len(regs))]

# dict to store the distances
dist_dict = {}

# get the nodata val
nodata_val = rio.open(hf.BIOCLIM_INFILEPATHS[0]).nodata

# columns for pandas DataFrame
seas_dist_colm = []
clim_dist_colm = []
kopp_dist_colm = []
reg_colm = []
x1_colm = []
x2_colm = []
y1_colm = []
y2_colm = []

if run_MMRR:
    MMRR_res = {}
if plot_ts:
    reg_ct = 1

ppt_file =  [f for f in hf.BIOCLIM_INFILEPATHS if 'bio_12.tif' in f]
assert len(ppt_file) == 1
ppt_file = ppt_file[0]
ppt = rxr.open_rasterio(ppt_file, masked=True)[0]

tmp_file =  [f for f in hf.BIOCLIM_INFILEPATHS if 'bio_1.tif' in f]
assert len(tmp_file) == 1
tmp_file = tmp_file[0]
tmp = rxr.open_rasterio(tmp_file, masked=True)[0]


mean_lats = []
mean_ppts = []
mean_tmps = []
for reg_poly, reg_pts, reg_col, reg_name in zip(regs,
                                                regs_pts,
                                                reg_cols,
                                                reg_names):

    print('\ngetting distance matrices for region: %s\n' % reg_name)

    # get points as nx2 numpy array
    pts = np.concatenate([np.array(pt.coords) for pt in reg_pts], axis=0)

    # get and plot the region's time series
    if plot_ts:
        tss = hf.get_raster_info_points(hf.COEFFS_FILE, pts, 'ts',
                                    standardize=standardize_ts)
        tss = tss[np.sum(np.isnan(tss), axis=1)==0, :]
        ax_ts = fig_ts.add_subplot(4,4,reg_ct)
        ax_ts.set_title(reg_name, fontdict={'fontsize': 10})
        for ts in tss:
            ax_ts.plot(ts, linestyle='-', color=reg_col, alpha=0.05)
        ax.set_ylim(np.min(tss), np.max(tss))
        reg_ct+=1

        # get points' pairwise clim dists
    clim_dist = hf.calc_pw_clim_dist_mat(pts, nodata_val=nodata_val)

    # get points' pairwise ts dists
    seas_dist = hf.get_raster_info_points(hf.COEFFS_FILE, pts, 'ts_pdm',
                                          standardize=standardize_ts)

    assert clim_dist.shape[0] == seas_dist.shape[0] == pts.shape[0]

    # get koppen differences for all points
    # (where 0 = same, 1 = different)
    kopp_vals = np.array([kopp.sel(x=x, y=y, method='nearest').values for x
                                                                    ,y in pts])

    kopp_dist = np.int8(distance.squareform(distance.pdist(
                                            np.atleast_2d(kopp_vals).T))>0)

    assert kopp_dist.shape[0] == clim_dist.shape[0]

    # drop clim dists for points without ts dists, and vice versa
    not_missing = np.where(np.nansum(seas_dist, axis=0)>0)[0]
    seas_dist = seas_dist[:, not_missing][not_missing,:]
    clim_dist = clim_dist[:, not_missing][not_missing,:]
    kopp_dist = kopp_dist[:, not_missing][not_missing,:]
    pts = pts[not_missing, :]
    still_not_missing = np.where(np.nansum(clim_dist, axis=0)>0)[0]
    seas_dist = seas_dist[:, still_not_missing][still_not_missing,:]
    clim_dist = clim_dist[:, still_not_missing][still_not_missing,:]
    kopp_dist = kopp_dist[:, still_not_missing][still_not_missing,:]
    pts = pts[still_not_missing, :]

    print(('\n%i points remain after dropping locations without seasonality '
          'data\n' % seas_dist.shape[0]))

    # run MMRR, if requested
    if run_MMRR:
        g = pyproj.Geod(ellps='WGS84')
        # get pw geo dist matrix of points
        geo_dist = np.ones([pts.shape[0]]*2) * np.nan
        for i in range(geo_dist.shape[0]):
            for j in range(i, geo_dist.shape[1]):
                lon1, lat1 = pts[i,:]
                lon2, lat2 = pts[j,:]
                (az12, az21, dist) = g.inv(lon1, lat1, lon2, lat2)
                geo_dist[i, j] = dist
                geo_dist[j, i] = dist
        # run model
        res = MMRR(seas_dist, [clim_dist, geo_dist], ['clim_dist', 'geo_dist'])
        # print and store results
        for k, v in res.items():
            print('\n\t%s: %s' % (k, str(v)))
        MMRR_res[reg_name] = res

    # extract the lower triangular values and scatter them
    indices = np.tril_indices(seas_dist.shape[0])

    seas_dist_vals = seas_dist[indices]
    clim_dist_vals = clim_dist[indices]
    kopp_dist_vals = kopp_dist[indices]

    # scatter 
    ax.scatter(clim_dist_vals, seas_dist_vals, s=1,
               c=reg_col, alpha=0.2, label=reg_name)

    # add the convex hull
    hull_pts = np.array((clim_dist_vals, seas_dist_vals)).T
    hull = ConvexHull(hull_pts)
    for simplex in hull.simplices:
        ax.plot(hull_pts[simplex, 0], hull_pts[simplex, 1], color=reg_col)

    # store the dists
    dist_dict[reg_name] = {'clim': clim_dist_vals,
                           'seas': seas_dist_vals,
                           'kopp': kopp_dist_vals,
                          }

    # add to DataFrame columns
    seas_dist_colm.extend(seas_dist_vals)
    clim_dist_colm.extend(clim_dist_vals)
    kopp_dist_colm.extend(kopp_dist_vals)
    reg_colm.extend([reg_name]*len(seas_dist_vals))
    x1_colm.extend(pts[indices[0],0])
    x2_colm.extend(pts[indices[1],0])
    y1_colm.extend(pts[indices[0],1])
    y2_colm.extend(pts[indices[1],1])

    # store the region's mean latitude and mean ppt
    mean_lat = [*reg_poly.centroid.coords.xy[1]][0]
    mean_ppt = float(np.mean(ppt.rio.clip(reg_poly)).values)
    mean_tmp = float(np.mean(tmp.rio.clip(reg_poly)).values)
    mean_lats.append(mean_lat)
    mean_ppts.append(mean_ppt)
    mean_tmps.append(mean_tmp)


ax.legend(fontsize=16)

ax.set_xlabel('climate distance (Euclidean)', size=18)
ax.set_ylabel('seasonal distance (Euclidean)', size=18)
ax.tick_params(labelsize=16)
fig.show()

# save the time-series plot grid
if plot_ts:
    fig_ts.subplots_adjust(bottom=0.1, top=0.92, left=0.1, right=0.92)
    fig_ts.savefig('time_series_plot_grid%s.png' % file_suffix, dpi=700)


# contourplot
df = pd.DataFrame({'seas_dist': seas_dist_colm,
                   'clim_dist': clim_dist_colm,
                   'kopp_dist': kopp_dist_colm,
                   'reg': reg_colm,
                   'x1': x1_colm,
                   'x2': x2_colm,
                   'y1': y1_colm,
                   'y2': y2_colm,
                  })

dists = []
for i, row in df.iterrows():
    dist = hf.calc_euc_dist(row[['x1', 'y1']].values, row[['x2', 'y2']].values)
    dists.append(dist)
df['geo_dist'] = dists

# write results to disk
df.to_csv('clim_dist_results_multiregion%s.csv' % file_suffix, index=False)
if run_MMRR:
    MMRR_res_df = pd.DataFrame(MMRR_res).T.reset_index()
    MMRR_res_df.columns = ['region']+[*MMRR_res_df.columns][1:]
    # add the mean latitude and mean ppt columns
    MMRR_res_df['mean_lat'] = mean_lats
    MMRR_res_df['mean_ppt'] = mean_ppts
    MMRR_res_df['mean_tmp'] = mean_tmps
    MMRR_res_df.to_csv('clim_dist_MMRR_results_multiregion%s.csv' % file_suffix, index=False)


# boxplot of same- versus different-Köppen-climate seas distances, by region
if make_boxplot:
    ax_plot = fig_map.add_subplot(gs[2:8, 1])
    sns.violinplot(x='kopp_dist', hue='reg', y='seas_dist', palette=reg_cols, data=df,
                   orient='h', ax=ax_plot)
    fig_map.subplots_adjust(top=0.96, bottom=0.1, left=0.07, right=0.97)
    fig_map.show()
    fig_map.savefig('kopp_seas_dist_results_multiregion%s.png' % file_suffix, dpi=700)

# save the map
fig_map.subplots_adjust(bottom=0.03, top=0.99, left=0.03, right=0.99)
fig_map.savefig('asynch_cluster_map_multiregion%s.png' % file_suffix, dpi=700)


# plot clim_dist MMRR coeff as a function of mean lat and mean ppt
fig_mmrr = plt.figure(figsize=(8,8))
ax_lat = fig_mmrr.add_subplot(1,2,1)
ax_lat.plot(MMRR_res_df['clim_dist'],
            MMRR_res_df['mean_lat'],
            linewidth=2,
            color='black',
            alpha=0.9,
            )
for r,x,y in zip(MMRR_res_df['region'].values,
                 MMRR_res_df['clim_dist'].values,
                 MMRR_res_df['mean_lat'].values):
    ax_lat.text(x, y, r, size=16, alpha=0.5)
ax_lat.axhline(y=0, color='red', linestyle=':')
ax_lat.set_xlabel('MMRR climate distance coefficient', fontdict={'fontsize': 18})
ax_lat.set_ylabel('mean latitude (degrees)', fontdict={'fontsize': 18})
ax_lat.tick_params(labelsize=14)


# plot on biomes
whittaker = pd.read_csv('./clim_dist/whittaker_biomes.csv', sep=';')

whittaker['temp_c'] = whittaker['temp_c'].apply(lambda x:
                                            float(x.replace(',', '.')))
whittaker['precp_mm'] = whittaker['precp_cm'].apply(lambda x:
                                            float(x.replace(',', '.'))*10)
biomes = []
centroids = []
patches = []

for biome in whittaker['biome'].unique():
    subwhit = whittaker[whittaker.biome == biome].loc[:, ['temp_c', 'precp_mm']].values
    centroids.append(np.mean(subwhit, axis=0))
    poly = mplPolygon(subwhit, True)
    patches.append(poly)
    biomes.append(re.sub('/', '/\n', biome))

#colors = ['#80fffd', # tundra
#          '#2b422f', # boreal forest
#          '#ebe157', # temperate grassland/desert
#          '#ab864b', # woodland/shrubland
#          '#17a323', # temperate seasonal forest
#          '#13916b', # temperate rain forest
#          '#00a632', # tropical rain forest
#          '#c2d69a', # tropical seasonal forest/savanna
#          '#e3a107', # subtropical desert
#         ]
#colors = np.array(colors)

colors = 255 * np.linspace(0, 1, len(patches))
p = PatchCollection(patches, alpha=0.4, edgecolor='k', cmap='Pastel1')
p.set_array(colors)
ax2 = fig_mmrr.add_subplot(1,2,2)
divider = make_axes_locatable(ax2)
cax2 = divider.append_axes('right', size='5%', pad=0.1)
ax2.add_collection(p)

for b,c in zip(biomes, centroids):
    ax2.text(c[0], c[1], b)

for r,x,y in zip(MMRR_res_df['region'].values,
                 MMRR_res_df['mean_tmp'].values,
                 MMRR_res_df['mean_ppt'].values):
    ax2.text(x, y, r, size=16, alpha=0.5)

scat = ax2.scatter(MMRR_res_df['mean_tmp'], MMRR_res_df['mean_ppt'],
                   c=MMRR_res_df['clim_dist'],
                   cmap='plasma_r')

plt.colorbar(scat, cax=cax2)
ax2.set_xlabel('MAT ($^{\circ}C$)',
              fontdict={'fontsize': 18})
ax2.set_ylabel('MAP ($mm$)',
              fontdict={'fontsize': 18})
ax2.tick_params(labelsize=14)

plt.show()

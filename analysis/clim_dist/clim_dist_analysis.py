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
from sklearn import cluster
from shapely.ops import unary_union
from shapely.geometry import Polygon, Point, MultiPolygon
from scipy import stats
from scipy.spatial import ConvexHull
from scipy.spatial import distance
from collections import Counter as C

# local modules
import helper_fns as hf
from MMRR import MMRR


# CLUSTER PARAMS:
# dbscan params
dbscan_eps = 3
dbscan_min_samples_frac_all = 0.5# alpha hull param
alpha = 1.25

# ANALYSIS PARAMS:
# standardize the seasonal time series before calculating their distances?
standardize_ts = True
# run MMRR tests?
run_MMRR = True


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

fig_map = plt.figure(figsize=(16,8))
ax_map = fig_map.add_subplot(111)
cntry.to_crs(4326).plot(color='none', edgecolor='k', linewidth=0.25, ax=ax_map)
ax_map.scatter(coords[:, 0], coords[:, 1], c=db.labels_, cmap='tab20', s=3)
ax_map.scatter(coords[:, 0], coords[:, 1], c=db.labels_>=0, cmap='Greys_r',
               alpha=db.labels_<0, s=3)

# get list of regions (i.e., cluster polygons)
regs = []
for clust_i in np.unique(db.labels_):
    # ignore noisy samples
    if clust_i >= 0:
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

# save the map
fig_map.subplots_adjust(bottom=0.03, top=0.99, left=0.03, right=0.99)
fig_map.savefig('asynch_cluster_map_multiregion.png', dpi=700)



####################################
# PART II: DRAW PTS AND RUN ANALYSIS
####################################

n_pts = 1000  # Enter in the number greater than random points you need

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

for reg_poly, reg_pts, reg_col, reg_name in zip(regs,
                                                regs_pts,
                                                reg_cols,
                                                reg_names):

    print('\ngetting distance matrices for region: %s\n' % reg_name)

    # get points as nx2 numpy array
    pts = np.concatenate([np.array(pt.coords) for pt in reg_pts], axis=0)

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

ax.legend(fontsize=16)

ax.set_xlabel('climate distance (Euclidean)', size=18)
ax.set_ylabel('seasonal distance (Euclidean)', size=18)
ax.tick_params(labelsize=16)
fig.show()


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
df.to_csv('clim_dist_results_multiregion.csv', index=False)
if run_MMRR:
    MMRR_res_df = pd.DataFrame(MMRR_res).T.reset_index()
    MMRR_res_df.columns = ['region']+[*MMRR_res_df.columns][1:]
    MMRR_res_df.to_csv('clim_dist_MMRR_results_multiregion.csv', index=False)


# boxplot of same- versus different-Köppen-climate seas distances, by region
fig_box = plt.figure(figsize=(16,8))
ax= fig_box.add_subplot(111)
sns.boxenplot(x='kopp_dist', hue='reg', y='seas_dist', palette=reg_cols, data=df,
            ax=ax)
fig_box.subplots_adjust(top=0.96, bottom=0.1, left=0.07, right=0.97)
fig_box.show()
fig_box.savefig('kopp_seas_dist_results_multiregion.png', dpi=700)

# and t-tests
#for lat_region in ['temperate', 'tropical']:
#    same = df[(df['kopp_dist']==0)&(df['reg']=='high asynch: %s' % lat_region)]['seas_dist']
#    diff = df[(df['kopp_dist']==1)&(df['reg']=='high asynch: %s' % lat_region)]['seas_dist']
#    mod = stats.ttest_ind(same, diff)
#    stat = mod.statistic
#    p = mod.pvalue
#    print('\n\nLATITUDINAL REGION: %s\n\n' % lat_region)
#    print('\n\nNUM SAME-KÖPPEN POINTS: %i\n\n' % len(same))
#    print('\n\nNUM DIFF-KÖPPEN POINTS: %i\n\n' % len(diff))
#    print('\n\nTEST RESULTS:\n\t%s\n\n' % str(mod))
#    print('-'*80+'\n\n')
#
#    # do permutations
#    kopp_dist_vals = np.array(deepcopy(kopp_dist_colm))
#    stat_list = []
#    p_list = []
#
#    for i in range(50):
#        np.random.shuffle(kopp_dist_vals)
#        same = df[(kopp_dist_vals==0)&(df['reg']=='high asynch: %s' % lat_region)]['seas_dist']
#        diff = df[(kopp_dist_vals==1)&(df['reg']=='high asynch: %s' % lat_region)]['seas_dist']
#        mod_perm = stats.ttest_ind(same, diff)
#        stat_list.append(mod_perm.statistic)
#        p_list.append(mod_perm.pvalue)
#
#    print('\n\nEMPIRICAL P-VALUE: %0.5e\n\n' % np.mean([np.abs(s)>=np.abs(stat)
#                                                        for s in stat_list]))


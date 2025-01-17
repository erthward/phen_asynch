#!/usr/bin/env python
# coding: utf-8

import geopandas as gpd
import numpy as np
import glob
import json
import time
import os
import sys
import re
import random
import pyproj
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
from copy import deepcopy
from scipy.spatial import KDTree
import statsmodels.api as sm
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
import statsmodels.api as sm
from scipy.spatial import ConvexHull
from scipy.spatial import distance
from collections import Counter as C
import itertools
import warnings
import h3

# local modules
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
from MMRR import MMRR
import phen_helper_fns as phf


print(f"\n\nSETTING PARAMETERS AND LOADING DATA...\n\n")
# BEHAVIORAL PARAMS:
# save CSVs of results?
save_all_results = True
# print everything?
verbose = True

# percentile above which to keep asynchrony pixels for clustering
top_asynch_pctile = 95

# CLUSTER PARAMS:
# dbscan params
dbscan_eps_vals = [2, 3.5, 5]
dbscan_minsamp_vals = [0.6, 0.45, 0.3]
# alpha hull param
alpha_vals = [0.25, 0.75, 1.25]

# minimum number of clusters to use for an analysis?
# (arbitrary, but nonetheless reasonable)
min_n_clusts = 7

# ANALYSIS PARAMS:
# asynch neighborhood radius (in km)
neigh_rad = 100

# a value that will wind up <= the number of random points
# used in the MMRR models (because some will fall in masked pixels)
n_pts = 1000

# standardize the seasonal time series before calculating their distances?
standardize_ts = True

# number of MMRR permutations
MMRR_nperm=499

# plot-formatting params
axlabel_fontdict = {'fontsize': 23}
ticklabel_size = 16
partlabel_size = 36



#####################################
# PART I: CLUSTER HIGH-ASYNCH REGIONS
#####################################

# load asynch data (band 0 in the asynch file
asynch = rxr.open_rasterio(phf.ASYNCH_FILES[neigh_rad])[0]

# load climate data
ppt_file =  [f for f in phf.BIOCLIM_INFILEPATHS if 'bio_12.tif' in f]
assert len(ppt_file) == 1
ppt_file = ppt_file[0]
ppt = rxr.open_rasterio(ppt_file, masked=True)[0]

tmp_file =  [f for f in phf.BIOCLIM_INFILEPATHS if 'bio_1.tif' in f]
assert len(tmp_file) == 1
tmp_file = tmp_file[0]
tmp = rxr.open_rasterio(tmp_file, masked=True)[0]

# set the parameterizations to loop over
loop_vals = itertools.product(dbscan_eps_vals,
                              dbscan_minsamp_vals,
                              alpha_vals)

# randomize order of the loop vals
# (so that I can rerun to try different parameterizations
#  and will eventually hit them all)
indices = [*range(len([*deepcopy(loop_vals)]))]
np.random.shuffle(indices)
loop_vals = [[*deepcopy(loop_vals)][i] for i in indices]

# get rid of loop vals that already exist in
# Shapefile containing partial output, if it exists
if os.path.isfile('clim_dep_all_MMRR_results_%ikmrad.shp' % neigh_rad):
    partial_results = gpd.read_file('clim_dep_all_MMRR_results_%ikmrad.shp' %
                                    neigh_rad)
    # map fiona's laundered column names back to original names
    partial_results = partial_results.rename(columns={
                                    'Intercept(': 'Intercept(t)',
                                    'clim_dist(': 'clim_dist(t)',
                                    'geo_dist(t': 'geo_dist(t)',
                                    'Intercep_1': 'Intercept(p)',
                                    'clim_dis_1': 'clim_dist(p)',
                                    'geo_dist(p': 'geo_dist(p)',
                                    'F-statisti': 'F-statistic',
                                    'mean_clim_': 'mean_clim_dist',
                                    'n_polys_dr': 'n_polys_dropped',
                                    'n_polys_al': 'n_polys_alpha_adjusted',
                                   })
    full_loop_vals = deepcopy(loop_vals)
    needed_loop_vals = []
    for dbscan_eps, dbscan_minsamp, alpha in deepcopy(loop_vals):
        subdf = partial_results[(partial_results['eps'] == dbscan_eps) &
                                (partial_results['minsamp'] == dbscan_minsamp) &
                                (partial_results['alpha'] == alpha)]
        if len(subdf) == 0:
            needed_loop_vals.append((dbscan_eps, dbscan_minsamp, alpha))
    print('\n\nOLD LENGTH LOOP VALS: %i' % len([*deepcopy(loop_vals)]))
    remaining_n_loop_vals = len([*deepcopy(needed_loop_vals)])
    print('\n\nNEW LENGTH LOOP VALS: %i\n\n' % remaining_n_loop_vals)
    loop_vals = needed_loop_vals
else:
    partial_results = None
    remaining_n_loop_vals = len([*deepcopy(loop_vals)])
    # start a new file to track all alpha-param adjustments for alpha hulls
    # that fail to form at the default alpha value tried
    # NOTE: this appears to happen when the dbscan algo throws the following
    # error: "WARNING:root:Singular matrix. Likely caused by all points lying
    # in an N-1 space."
    with open('clim_dep_alpha_adjust_log.csv', 'w') as f:
        f.write('loop_num,eps,minsamp,alpha,clust_i,clust_mean_x,clust_mean_y,alpha_adustment')
    print('\n\nNO EXISTING RESULTS DETECTED\n\n')

# create data strucutres to save all results for final analysis
all_loop_coords = {}
all_loop_polys = {}
all_loop_rand_pts = {}
all_loop_MMRR_res = {}

# add previously saved partial results to the all_loop_MMRR_res dict,
# if there were any
if partial_results is not None:
    prev_param_combos = partial_results[['alpha', 'minsamp', 'eps']].values
    n_prev_param_combos = len(set([tuple(row) for row in prev_param_combos]))
    print('\n\n%i PREVIOUS PARAM COMBOS ALREADY RUN' % n_prev_param_combos)
    # NOTE: MAKE LOOP NUM -1 FOR ALL PREV RESULTS, EVEN IF MULTIPLE LOOPS'
    #       RESULTS, SINCE THE LOOP NUM IS NOT REALLY USEFUL POST HOC ANYHOW
    all_loop_MMRR_res[-1] = partial_results
else:
    n_prev_param_combos = 0

# mask everything below Nth percentile
pctile = np.nanpercentile(asynch, top_asynch_pctile)
maxes = asynch.where(asynch>=pctile, np.nan).clip(min=1, max=1)
# extract as points
X, Y = np.meshgrid(maxes.x.values, maxes.y.values)
coords = np.array([*zip(X.ravel(), Y.ravel())])
coords = coords[np.where(pd.notnull(maxes.values.ravel()))[0]]

# loop analysis over param vals
loop_ct = 0
if remaining_n_loop_vals > 0:
    for (dbscan_eps, dbscan_minsamp, alpha) in loop_vals:
        start_time = time.time()
        n_polys_dropped = 0
        n_polys_alpha_adjusted = 0
        print((('#'*80+'\n')*4).rstrip())
        pct_done = 100 * (loop_ct+n_prev_param_combos+1)/(len(dbscan_eps_vals) *
                                                        len(dbscan_minsamp_vals) *
                                                        len(alpha_vals))
        print((f'LOOP {loop_ct} '
               f'({np.round(pct_done, 1)}%):'))
        print('#'*16 + '\n\n')
        print('\tdbscan_eps: %0.3f\n\n' % dbscan_eps)
        print('\tdbscan_minsamp: %0.3f\n\n' % dbscan_minsamp)
        print('\talpha: %0.3f\n\n' % alpha)

        print(f"\n\nCLUSTERING HIGH-ASYNCHRONY REGIONS...\n\n")
        # run clustering algo
        # NOTE: DBSCAN SEEMED LIKE THE BEST FIT OFF THE BAT,
        #       BASED ON DEPICTION HERE:
        #           https://scikit-learn.org/stable/modules/clustering.html#overview-of-clustering-methods
        #       AND DESCRIPTION HERE:
        #           https://scikit-learn.org/stable/modules/clustering.html#dbscan
        db = cluster.DBSCAN(eps=dbscan_eps,
                min_samples=dbscan_minsamp * maxes.shape[0]).fit(coords)
        # NOTE: ADD 1 TO CLUSTER LABELS, SO THAT FIRST CLUSTER == 1, NOT 0,
        #       AND THIS WAY THE 'NOISY' SAMPLES THAT DBSCAN DIDN'T INCORPORATE
        #       INTO ANY CLUSTERS WILL BE LABELED 0, NOT -1
        db.labels_ += 1

        # get list of regions (i.e., cluster polygons)
        polys = []
        for clust_i in np.unique(db.labels_):
            # ignore the 'noisy' samples (i.e., those labeled 0)
            if clust_i > 0:
                if verbose:
                    print('\n\n\tgetting polygon and mapping for cluster %i...\n\n' % clust_i)
                # indices of this cluster's samples
                inds = np.where(db.labels_==clust_i)[0]
                # inds = [*set([*np.where(db.labels_==clust_i)[0]]).intersection(
                         #set(db.core_sample_indices_))]
                # coords of this cluster's samples
                clust_coords = coords[inds,:]
                # get alpha-hull coords
                # NOTE: looping the alphashape function and adding 0.05 to alpha each
                # time because I ran into some weird issues with the function where for
                # certain point sets (e.g., South African cape one specific time) it
                # returned an empty GeometryCollection rather than a Poly or
                # MultiPoly...
                alpha_hull = None
                ct = 0
                while not (isinstance(alpha_hull, Polygon) or isinstance(alpha_hull,
                                                                         MultiPolygon)):
                    alpha_hull = alphashape.alphashape(clust_coords,
                                                       alpha=alpha+(0.05*ct))
                    ct+=1
                if ct > 1:
                    n_polys_alpha_adjusted += 1
                    print('\n\n\tNOTE: ADDED 0.05 * %i TO ALPHA!\n\n' % (ct-1))
                    with open('clim_dep_alpha_adjust_log.csv', 'a') as f:
                        f.write((f"\n{loop_ct},{dbscan_eps},{dbscan_minsamp},{alpha},{clust_i},"
                                 f"{np.mean(clust_coords[:,0])},{np.mean(clust_coords[:,1])},"
                                 f"{(ct-1)*0.05}")
                               )
                # coerce to MultiPolygon, if necessary
                if isinstance(alpha_hull, MultiPolygon):
                    pass
                else:
                    alpha_hull = MultiPolygon([alpha_hull])
                # save the MultiPolygon
                polys.append(alpha_hull)

        # coerce all to MultiPolygons (to be able to use for raster clipping
        polys = [poly if isinstance(poly,
                    MultiPolygon) else MultiPolygon([poly]) for poly in polys]

        # drop all sliver polygons
        # (i.e., polygons with projected areas < equivalent of ~5x5 pixels in
        # our raster, or < (5500m *5)**2 m^2 area),
        # as these are the ones that cause the draw_random_points fn's
        # while loop to run forever, since points never fall inside,
        # and because they'd also be meaningless to analyze anyhow
        len_b4 = len(polys)
        min_area = (5500*5)**2
        areas = gpd.GeoDataFrame(pd.DataFrame({'a':[*range(len(polys))],
                                               'polys': polys}),
                                 geometry='polys',
                                 crs=4326,
                                ).to_crs(8857).area
        polys = [p for p, a in zip(polys, areas) if a >= min_area]
        len_af = len(polys)
        n_polys_dropped += len_b4-len_af
        if verbose and len_af != len_b4:
            print((f"\n\tNOTE: dropped {n_polys_dropped} tiny polygon"
                   f"{'s' * ((len_b4-len_af)>1)} from analysis\n"))


####################################
# PART II: DRAW PTS AND RUN ANALYSIS
####################################

        print(f"\n\nRUNNING REGIONAL ANALYSES...\n\n")
        print(f"\n\n drawing random points...\n\n")
        # NOTE: drawing considerably more than desired, 
        #       so that we get closer to the target
        #       even after attrition because of masked LSP pixels
        polys_pts = [phf.draw_random_points_within_polygon(
                            3*n_pts, poly, verbose) for poly in polys]

        # list of region names
        reg_names = [str(i+1) for i in range(len(polys))]

        # create colors for the regions
        reg_cols = [plt.cm.tab20(v/len(polys)) for v in range(len(polys))]

        # get the nodata val
        nodata_val = rio.open(phf.BIOCLIM_INFILEPATHS[0]).nodata

        # data structure for the MMRR results
        MMRR_res = {}

        mean_lats = []
        mean_ppts = []
        mean_tmps = []
        mean_clim_dists = []

        # run an MMRR model for each delineated region
        for reg_poly, reg_pts, reg_col, reg_name in zip(polys,
                                                        polys_pts,
                                                        reg_cols,
                                                        reg_names):
            if verbose:
                print('\ngetting distance matrices for region: %s\n' % reg_name)

            # get points as nx2 numpy array
            pts = reg_pts['pts'].get_coordinates().values
            # get points' pairwise clim dists
            clim_dist = phf.calc_pw_clim_dist_mat(pts, nodata_val=nodata_val)

            # get points' pairwise ts dists
            seas_dist = phf.get_raster_info_points(phf.COEFFS_STRICT_FILE, pts, 'ts_pdm',
                                                   standardize=standardize_ts)

            assert clim_dist.shape[0] == seas_dist.shape[0] == pts.shape[0]

            # drop clim dists for points without ts dists, and vice versa
            not_missing = np.where(np.nansum(seas_dist, axis=0)>0)[0]
            seas_dist = seas_dist[:, not_missing][not_missing,:]
            clim_dist = clim_dist[:, not_missing][not_missing,:]
            pts = pts[not_missing, :]
            still_not_missing = np.where(np.nansum(clim_dist, axis=0)>0)[0]
            seas_dist = seas_dist[:, still_not_missing][still_not_missing,:]
            clim_dist = clim_dist[:, still_not_missing][still_not_missing,:]
            pts = pts[still_not_missing, :]
            # only keep n_pts number of points (in case more remain)
            assert (seas_dist.shape[0] == seas_dist.shape[1] ==
                    clim_dist.shape[0] == clim_dist.shape[1] == pts.shape[0])
            if pts.shape[0] > n_pts:
                seas_dist = seas_dist[:n_pts, :n_pts]
                clim_dist = clim_dist[:n_pts, :n_pts]
                pts = pts[:n_pts, :]
            if verbose:
                print(('\n\t%i points remain after dropping locations '
                       'without seasonality data\n') % seas_dist.shape[0])
                print(f"\n\n\trunning MMRR...\n\n")
            if seas_dist.shape[0]>=2:

                g = pyproj.Geod(ellps='WGS84')
                # get pw geo dist matrix of points
                geo_dist= np.ones([pts.shape[0]]*2) * np.nan
                for i in range(geo_dist.shape[0]):
                    for j in range(i, geo_dist.shape[1]):
                        lon1, lat1 = pts[i,:]
                        lon2, lat2 = pts[j,:]
                        (az12, az21, dist) = g.inv(lon1, lat1, lon2, lat2)
                        geo_dist[i, j] = dist
                        geo_dist[j, i] = dist

                # run model
                res = MMRR(Y=seas_dist,
                           X=[clim_dist, geo_dist],
                           Xnames=['clim_dist', 'geo_dist'],
                           # NOTE: MMRR will standardize lower-triangular
                           #       values, and thus provide coeff values as
                           #       beta-coefficients
                           standardize=True,
                           intercept=True,
                           nperm=MMRR_nperm,
                          )

                # store results
                if verbose:
                    print(f"\n\n\tstoring results...\n\n")
                MMRR_res[reg_name] = res

                # extract the lower triangular values
                indices = np.tril_indices(seas_dist.shape[0])

                seas_dist_vals = seas_dist[indices]
                clim_dist_vals = clim_dist[indices]
                geo_dist_vals = geo_dist[indices]

                # store the polygon's mean latitude and mean climate vals
                mean_lat = [*reg_poly.centroid.coords.xy[1]][0]
                # NOTE: ignore warning about deprecated iteration of multi-part
                #       shapely geometries, for now
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    mean_ppt = float(np.mean(ppt.rio.clip([reg_poly])).values)
                    mean_tmp = float(np.mean(tmp.rio.clip([reg_poly])).values)
                mean_clim_dist = np.mean(clim_dist_vals)
                mean_lats.append(mean_lat)
                mean_ppts.append(mean_ppt)
                mean_tmps.append(mean_tmp)
                mean_clim_dists.append(mean_clim_dist)

        MMRR_res_df = pd.DataFrame(MMRR_res).T.reset_index()
        MMRR_res_df.columns = ['region']+[*MMRR_res_df.columns][1:]

        # add the mean latitude and climate columns
        MMRR_res_df['mean_lat'] = mean_lats
        MMRR_res_df['mean_ppt'] = mean_ppts
        MMRR_res_df['mean_tmp'] = mean_tmps
        MMRR_res_df['mean_clim_dist'] = mean_clim_dists

        # add the loop_ct counter and param vals as columns, so that all
        # saved results dfs can easily be concatenated later
        MMRR_res_df['loop_num'] = loop_ct
        MMRR_res_df['eps'] = dbscan_eps
        MMRR_res_df['minsamp'] = dbscan_minsamp
        MMRR_res_df['alpha'] = alpha

        # add some performance metrics
        MMRR_res_df['n_polys_dropped'] = n_polys_dropped
        MMRR_res_df['n_polys_alpha_adjusted'] = n_polys_alpha_adjusted
        end_time = time.time()
        runtime = end_time-start_time
        MMRR_res_df['runtime'] = runtime

        # add the polygons as a column
        MMRR_res_df['geometry'] = polys

        # convert to a GeoDataFrame
        MMRR_res_gdf = gpd.GeoDataFrame(MMRR_res_df, geometry='geometry', crs=4326)

        # check that columns line up with partial_results columns,
        # if partial results are being used
        if partial_results is not None:
            assert np.all(MMRR_res_gdf.columns == partial_results.columns)

        # store all other results for this loop
        all_loop_coords[loop_ct] = coords
        all_loop_polys[loop_ct] = polys
        all_loop_rand_pts[loop_ct] = polys_pts
        all_loop_MMRR_res[loop_ct] = MMRR_res_gdf

        # concatenate all MMRR results' dataframes
        # (do this each loop so that if the script breaks partway then I can
        # pick up where I left off)
        all_loop_MMRR_res_gdf = pd.concat([df.to_crs(4326) for df in all_loop_MMRR_res.values()])
        # check columns all got lined up correctly
        assert np.all(all_loop_MMRR_res_gdf.columns == MMRR_res_gdf.columns)
        if save_all_results:
            all_loop_MMRR_res_gdf.to_file(
                'clim_dep_all_MMRR_results_%ikmrad.shp' % neigh_rad,
                                          index=False)

        # increment loop count
        loop_ct += 1

# if all results already run (i.e., if remaining_n_loop_vals == 0),
# then just prep data for analysis
else:
    all_loop_MMRR_res_gdf = partial_results



####################
# PART III: ANALYSIS
####################
print(f"\n\nRUNNING OVERALL ANALYSES...\n\n")
fig = plt.figure(figsize=(23,7))
gs = fig.add_gridspec(nrows=100, ncols=300)

# run single, overarching OLS model of clim_dist coeff as a function of mean
# cluster latitude
x = sm.add_constant(np.abs(all_loop_MMRR_res_gdf['mean_lat']))
y = all_loop_MMRR_res_gdf['clim_dist']
mod = sm.OLS(endog=y, exog=x).fit()
coeffs = mod.params
pvalues = mod.pvalues
fx_size = coeffs.mean_lat

# because results from multiple MMRR models are not independent (since clusters
# overlap in space), perform a permutation test of the 
null_fx = []
for _ in range(1000):
    null_x = sm.add_constant(np.abs(np.random.choice(all_loop_MMRR_res_gdf['mean_lat'],
                                     size=len(all_loop_MMRR_res_gdf),
                                     replace=False)))
    null_mod = sm.OLS(endog=y, exog=null_x).fit()
    null_coeffs = null_mod.params
    null_fx_size = null_coeffs.x1
    null_fx.append(null_fx_size)

# display empirical p-value (i.e., percent of null results >= real result)
emp_pval = np.mean(np.array(null_fx) >= fx_size)
# NOTE: fix empirical p-value at 0.001 for display, if == 0
if emp_pval == 0:
    emp_pval = 0.001
emp_pval_str = '\n\nEMPIRICAL P-VALUE: %0.3f\n' % emp_pval
print(emp_pval_str)


#........................
# PART C: empirical p-val
#........................

# plot clim_dist effect size vs. empirical null dist
ax = fig.add_subplot(gs[80:, 215:])
ax.hist(null_fx, bins=100, alpha=0.8, color='black')
ax.axvline(fx_size, *ax.get_ylim(), color='red', linewidth=2.5)
ax.set_xlabel('slope', fontdict=axlabel_fontdict)
ax.set_ylabel('freq.', fontdict=axlabel_fontdict)
xticks = np.linspace(ax.get_xticks()[0], ax.get_xticks()[-1], 5)
ax.set_xticks(xticks, ['%0.3f' % t for t in xticks])
ax.set_yticks(())
ax.tick_params(labelsize=ticklabel_size)
#ax.text(0.2*ax.get_xlim()[1],
#        0.75*ax.get_ylim()[1],
#        'observed = %0.5f' % fx_size,
#        color='red',
#        size=16,
#       )
ax.text(0.35*ax.get_xlim()[1],
        0.75*ax.get_ylim()[1],
        '$P_{permut} \ll %s' % f'{np.round(emp_pval, 3)}$',
        color='black',
        size=16,
       )

# add label for part C
ax.text(ax.get_xlim()[0] - (0.2*(ax.get_xlim()[1]-ax.get_xlim()[0])),
        1.05*ax.get_ylim()[1],
        'C.', size=partlabel_size, weight='bold')


#.....................
# PART B: scatter plot
#.....................

# plot clim_dist effect size vs |mean lat|
ylabel = 'effect size of mean cluster clim_dist\non mean cluster seas_dist'
ax = fig.add_subplot(gs[:60, 215:])
ax.scatter(np.abs(all_loop_MMRR_res_gdf['mean_lat']),
            all_loop_MMRR_res_gdf['clim_dist'],
            c='black', s=35)
xlim = np.array(ax.get_xlim())
ax.plot(xlim,
        coeffs.const + coeffs.mean_lat * xlim,
        linewidth=2.5,
        alpha=0.7,
        color='red',
       )
ax.set_xlabel('mean absolute latitude', fontdict=axlabel_fontdict)
ax.set_ylabel('mean $β_c$', fontdict=axlabel_fontdict)
ax.tick_params(labelsize=ticklabel_size)

# add label for part B
ax.text(ax.get_xlim()[0] - (0.2*(ax.get_xlim()[1]-ax.get_xlim()[0])),
        1.05*ax.get_ylim()[1],
        'B.', size=partlabel_size, weight='bold')
# add slope and p-value
slope = coeffs.mean_lat
pval = pvalues.mean_lat
# NOTE: fix p-value at 0.001 for display, if == 0
if np.round(pval, 3) == 0:
    pval = 0.001
ax.text(-1.5,
        0.9*ax.get_ylim()[1],
        'slope = %0.5f' % slope,
        color='red',
        size=16,
       )
ax.text(-1.5,
        0.78*ax.get_ylim()[1],
        f'$P \ll {np.round(pval, 3)}$',
        color='black',
        size=16,
       )



#........................
# PART A: summary hex map
#........................

# NOTE: based on code taken from:
    # https://towardsdatascience.com/uber-h3-for-data-analysis-with-python-1e54acdcc908
ax = fig.add_subplot(gs[:, :190])

# add bottom axes for a colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes('bottom', size='7%', pad=0.2)

# make dataframe to hold h3-converted data
h3_df = pd.DataFrame([], columns=['row_id', 'h3_id',
                                  'h3_geo_boundary', 'h3_centroid'])

# loop over results rows and convert to H3 hexes
for i, row in all_loop_MMRR_res_gdf.iterrows():
    p = row['geometry']
    if isinstance(p, MultiPolygon):
        ps = [*p.geoms]
    else:
        ps = [p]
    for poly in ps:
        poly_json = gpd.GeoSeries([poly]).__geo_interface__
        poly_json = poly_json['features'][0]['geometry']
        h3_hexes = h3.polyfill_geojson(poly_json, 3)
        for h3_hex in h3_hexes:
            h3_geo_boundary = Polygon(
                h3.h3_to_geo_boundary(h3_hex,geo_json=True))
            h3_centroid = h3.h3_to_geo(h3_hex)
            h3_df.loc[len(h3_df)] = [i, h3_hex, h3_geo_boundary, h3_centroid]

# coerce to GeoDataFrame
geoms = [Polygon(row['h3_geo_boundary']) for i,
                                            row in h3_df.iterrows()]
h3_df['geometry'] = geoms
h3_gdf = gpd.GeoDataFrame(h3_df, geometry='geometry', crs=4326)

# deduplicate hexes
h3_gdf = h3_gdf.drop_duplicates(subset='geometry')

# summarize results within hexes
mean_clim_dist_vals = []
region_counts = []
for i, row in h3_gdf.iterrows():
    mean_clim_dist_val = float(all_loop_MMRR_res_gdf[
        all_loop_MMRR_res_gdf.intersects(row['geometry'])].loc[:,
                                                               'clim_dist'].mean())
    region_count = float(all_loop_MMRR_res_gdf[
        all_loop_MMRR_res_gdf.intersects(row['geometry'])].loc[:,
                                                               'clim_dist'].count())
    mean_clim_dist_vals.append(mean_clim_dist_val)
    region_counts.append(region_count)
assert len(mean_clim_dist_vals) == len(h3_gdf)
h3_gdf['clim_dist_mean'] = mean_clim_dist_vals
h3_gdf['scaled_counts'] = np.array(region_counts)/np.max(region_counts)

# transform to equal-area projection and plot
hex_subplot = h3_gdf.to_crs(8857).plot('clim_dist_mean',
                                       cmap='magma',
                                       alpha=0.8,
                                       zorder=2,
                                       ax=ax,
                                       edgecolor='white',
                                       linewidth=0.2,
                                      )
phf.plot_juris_bounds(ax,
                      lev1_alpha=0.7,
                      lev1_zorder=3,
                      lev1_linewidth=0.05,
                      lev0_alpha=1,
                      lev0_zorder=4,
                      lev0_linewidth=0.15,
                      strip_axes=True,
                     )

# add label for part A
ax.text(1.1*ax.get_xlim()[0],
        1.06*ax.get_ylim()[1],
        'A.', size=partlabel_size, weight='bold')

# add the colorbar
# NOTE: need to create a custom ScalarMappable to feed into the colormap call
scalmap = plt.cm.ScalarMappable(cmap='magma',
                        norm=plt.Normalize(vmin=np.min(h3_gdf.clim_dist_mean),
                                           vmax=np.max(h3_gdf.clim_dist_mean)))
plt.colorbar(scalmap, cax=cax, orientation='horizontal')
xticks = np.linspace(np.min(h3_gdf.clim_dist_mean),
                     np.max(h3_gdf.clim_dist_mean), 5)

cax.set_xlabel('mean $β_c$', fontdict=axlabel_fontdict)


cax.set_xticks(xticks, ['%0.2f' % t for t in xticks], size=ticklabel_size)
cax.set_ylabel('')
cax.set_yticks(())



# adjust subplots and write to disk
fig.subplots_adjust(top=0.92,
                    bottom=0.125,
                    right=0.98,
                    left=0.026,
                    wspace=0.5,
                    hspace=0.3,
                   )

fig.savefig(os.path.join(phf.FIGS_DIR, 'FIG_isoclim_asynch.png'), dpi=700)


#!/usr/bin/env python
# coding: utf-8

"""
# TODO:

    - include Whittaker plot?

    - attend remaining TODO's

    - need to record places where alpha is increased, or get rid of that, or is
         that okay as part of clustering?

    - need to weed out clim_dist coeff values when their p-values are >0.05, or
         is that part of the point? (I looked and they're mostly far-north places
         that will go away with proper masking anyhow...)
"""

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
import statsmodels.api as sm
from scipy.spatial import ConvexHull
from scipy.spatial import distance
from collections import Counter as C
import itertools
import warnings

# local modules
from MMRR import MMRR
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf

# BEHAVIORAL PARAMS:
# plot Whittaker biomes?
plot_whittaker = True
# save CSVs of results?
save_all_results = True
# print everything?
verbose = False

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
# asynch neighborhood radius (in m)
neigh_rad = 100_000

# a number <= number of random points desired
n_pts = 1000

# standardize the seasonal time series before calculating their distances?
standardize_ts = True

# number of MMRR permutations
MMRR_nperm=499



#####################################
# PART I: CLUSTER HIGH-ASYNCH REGIONS
#####################################

# load country boundaries
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                       'NewWorldFile_2020.shp')).to_crs(8857)

# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR,f)).to_crs(8857))
subnational = pd.concat(subnational)

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

# determine if in interactive mode,
# then set accordingly the parameterizations to loop over
interactive = hasattr(sys, 'ps1')
if interactive:
    loop_vals = itertools.product(dbscan_eps_vals,
                                  dbscan_minsamp_vals,
                                  alpha_vals)
else:
    loop_vals = [[float(arg.strip()) for arg in sys.argv[1:]]]

# randomize order of the loop vals
# (so that I can rerun to try different parameterizations
#  and will eventually hit them all)
indices = [*range(len([*deepcopy(loop_vals)]))]
np.random.shuffle(indices)
loop_vals = [[*deepcopy(loop_vals)][i] for i in indices]

# get rid of loop vals that already exist in CSV of partial output,
# if it exists
if os.path.isfile('clim_dist_all_MMRR_results_%imrad.shp' % neigh_rad):
    partial_results = gpd.read_file('clim_dist_all_MMRR_results_%imrad.shp' %
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
    print('\n\nNO EXISTING RESULTS DETECTED\n\n')

# create data strucutres to save all results for final analysis
all_loop_coords = {}
all_loop_labels = {}
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

# loop analysis over param vals
loop_ct = 0
if remaining_n_loop_vals > 0:
    for (dbscan_eps, dbscan_minsamp, alpha) in loop_vals:
        print(('#'*80+'\n')*4)
        print(f'\n\nLOOP {loop_ct}:\n\n')
        print('\tdbscan_eps: %0.3f\n\n' % dbscan_eps)
        print('\tdbscan_minsamp: %0.3f\n\n' % dbscan_minsamp)
        print('\talpha: %0.3f\n\n' % alpha)

        # mask everything below Nth percentile
        pctile = np.nanpercentile(asynch, 95)
        maxes = asynch.where(asynch>=pctile, np.nan).clip(min=1, max=1)
        # extract as points
        X, Y = np.meshgrid(maxes.x.values, maxes.y.values)
        coords = np.array([*zip(X.ravel(), Y.ravel())])
        coords = coords[np.where(pd.notnull(maxes.values.ravel()))[0]]

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
                    print('\n\nGETTING POLYGON AND MAPPING FOR CLUSTER %i...\n\n' % clust_i)
                # indices of this cluster's samples
                inds = np.where(db.labels_==clust_i)[0]
                # inds = [*set([*np.where(db.labels_==clust_i)[0]]).intersection(
                         #set(db.core_sample_indices_))]
                # coords of this cluster's samples
                clust_coords = coords[inds,:]
                # get alpha-hull coords
                # TODO: DOESN'T THIS VIOLATE THE FACT THAT I'M CHOOSING SET
                #       ALPHA VALUES? FIGURE OUT WHAT TO DO HERE...
                # NOTE: looping the alphashape function and adding 0.25 to alpha each
                # time because I ran into some weird issues with the function where for
                # certain point sets (e.g., South African cape one specific time) it
                # returned an empty GeometryCollection rather than a Poly or
                # MultiPoly...
                alpha_hull = None
                ct = 0
                while not (isinstance(alpha_hull, Polygon) or isinstance(alpha_hull,
                                                                         MultiPolygon)):
                    alpha_hull = alphashape.alphashape(clust_coords,
                                                       alpha=alpha+(0.25*ct))
                    ct+=1
                if ct > 1:
                    print('\n\n\tNOTE: ADDED 0.25 * %i TO ALPHA!\n\n' % ct)
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



####################################
# PART II: DRAW PTS AND RUN ANALYSIS
####################################

        polys_pts = [phf.generate_random_points_in_polygon(
                                                n_pts, poly) for poly in polys]

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
        for reg_poly, reg_pts, reg_col, reg_name in zip(polys,
                                                        polys_pts,
                                                        reg_cols,
                                                        reg_names):
            if verbose:
                print('\ngetting distance matrices for region: %s\n' % reg_name)

            # get points as nx2 numpy array
            pts = np.concatenate([np.array(pt.coords) for pt in reg_pts], axis=0)
            # get points' pairwise clim dists
            clim_dist = phf.calc_pw_clim_dist_mat(pts, nodata_val=nodata_val)

            # get points' pairwise ts dists
            seas_dist = phf.get_raster_info_points(phf.COEFFS_FILE, pts, 'ts_pdm',
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
            if verbose:
                print(('\n%i points remain after dropping locations '
                       'without seasonality data\n') % seas_dist.shape[0])
            if seas_dist.shape[0]>=2:

                # run MMRR
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

                # standardize input vars (to make MMRR coeffs into standardized
                # beta coeffs, and thus comparable)
                seas_dist_stand = phf.standardize_array(seas_dist)
                clim_dist_stand = phf.standardize_array(clim_dist)
                geo_dist_stand = phf.standardize_array(geo_dist)
                assert np.allclose(np.std(seas_dist_stand), 1)
                assert np.allclose(np.std(clim_dist_stand), 1)
                assert np.allclose(np.std(geo_dist_stand), 1)

                # run model
                res = MMRR(seas_dist_stand,
                           [clim_dist_stand, geo_dist_stand],
                           ['clim_dist', 'geo_dist'],
                           MMRR_nperm)

                # store results
                MMRR_res[reg_name] = res

                # extract the lower triangular values and scatter them
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
                    mean_ppt = float(np.mean(ppt.rio.clip(reg_poly)).values)
                    mean_tmp = float(np.mean(tmp.rio.clip(reg_poly)).values)
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
        all_loop_labels[loop_ct] = db.labels_
        all_loop_polys[loop_ct] = polys
        all_loop_rand_pts[loop_ct] = polys_pts
        all_loop_MMRR_res[loop_ct] = MMRR_res_gdf

        # concatenate all MMRR results' dataframes
        # (do this each loop so that if the script breaks partway then I can
        # pick up where I left off)
        all_loop_MMRR_res_gdf = pd.concat([*all_loop_MMRR_res.values()])
        # check columns all got lined up correctly
        assert np.all(all_loop_MMRR_res_gdf.columns == MMRR_res_gdf.columns)
        if save_all_results:
            all_loop_MMRR_res_gdf.to_file(
                'clim_dist_all_MMRR_results_%imrad.shp' % neigh_rad,
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

# run single, overarching OLS model of clim_dist coeff as a function of mean
# cluster latitude
x = all_loop_MMRR_res_gdf['mean_lat']
y = np.abs(all_loop_MMRR_res_gdf['clim_dist'])
# TODO: WHETHER OR NOT TO INCLUDE INTERCEPT?
mod = sm.OLS(y, x).fit()
coeffs = mod.params
fx_size = coeffs[0]

# because results from multiple MMRR models are not independent (since clusters
# overlap in space), perform a permutation test of the 
null_fx = []
for _ in range(1000):
    null_x = np.random.choice(all_loop_MMRR_res_gdf['mean_lat'],
                              size=len(all_loop_MMRR_res_gdf),
                              replace=False)
    null_mod = sm.OLS(y, null_x).fit()
    null_coeffs = null_mod.params
    null_fx_size = null_coeffs[0]
    null_fx.append(null_fx_size)

# display empirical p-value (i.e., percent of null results >= real result)
emp_pval = '\n\nEMPIRICAL P-VALUE: %0.3f\n' % np.mean(np.array(null_fx) >= fx_size)
print(emp_pval)

# plot clim_dist effect size vs. empirical null dist
fig_null = plt.figure(figsize=(6,6))
ax = fig_null.add_subplot(111)
ax.set_xlabel('effect size of mean cluster latitude on MMRR clim_dist coeff')
ax.set_ylabel('freq')
plt.hist(null_fx, bins=100, alpha=0.7)
ax.axvline(fx_size, *ax.get_ylim(), color='red', linewidth=2)
ax.set_title(emp_pval)
fig_null.savefig('mean_cluster_lat_effect_size_null_dist_%imrad.png' % neigh_rad, dpi=700)

# plot clim_dist effect size vs |mean lat|
fig_scat = plt.figure(figsize=(6,12))
ylabels = {'clim_dist': ('effect size of mean cluster clim_dist\n'
                         'on mean cluster seas_dist'),
           'clim_dist(p)': 'p-value of clim_dist coefficient in MMRR',
           'mean_clim_dist': 'cluster mean climate distance',
          }
for i, y_col in enumerate(['clim_dist', 'clim_dist(p)', 'mean_clim_dist']):
    ax = fig_scat.add_subplot(311+i)
    if i == 2:
        ax.set_xlabel('|mean cluster latitude|')
    else:
        ax.set_xlabel('')
    ax.set_ylabel(ylabels[y_col])
    ax.scatter(np.abs(all_loop_MMRR_res_gdf['mean_lat']),
                all_loop_MMRR_res_gdf[y_col],
                c='black', s=35)
    sns.regplot(x=np.abs(all_loop_MMRR_res_gdf['mean_lat']),
                y=all_loop_MMRR_res_gdf[y_col],
                scatter=False,
                ax=ax,
               )

fig_scat.savefig('clim_dist_vs_lat_scat_%imrad.png' % neigh_rad, dpi=700)

fig_map = plt.figure(figsize=(20,10))
ax = fig_map.add_subplot(111)
subnational.to_crs(4326).plot(color='none', edgecolor='black', zorder=0, ax=ax,
                             alpha=0.6)
countries.to_crs(4326).plot(color='none', edgecolor='black', zorder=1, ax=ax)
gdf_for_plot = deepcopy(all_loop_MMRR_res_gdf)
gdf_for_plot.geometry = gdf_for_plot.boundary
gdf_for_plot.plot(column='clim_dist',
                  alpha=0.8,
                  #s=50,
                  #edgecolor='black',
                  linewidth=2,
                  zorder=2,
                  cmap='plasma',
                  ax=ax
                 )
fig_map.savefig('asynch_clusts_map_%imrad.png' % neigh_rad, dpi=700)

plt.show()

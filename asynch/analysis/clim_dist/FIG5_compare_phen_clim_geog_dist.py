#!/usr/bin/env python
# coding: utf-8

"""
# TODO:

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
import h3

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
# asynch neighborhood radius (in km)
neigh_rad = 100

# a number <= number of random points desired
n_pts = 1000

# standardize the seasonal time series before calculating their distances?
standardize_ts = True

# number of MMRR permutations
MMRR_nperm=499


# plot formatting params
axlabel_fontdict = {'fontsize': 25}
ticklabel_size = 16
partlabel_size = 36



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

# get rid of loop vals that already exist in CSV of partial output,
# if it exists
if os.path.isfile('clim_dist_all_MMRR_results_%ikmrad.shp' % neigh_rad):
    partial_results = gpd.read_file('clim_dist_all_MMRR_results_%ikmrad.shp' %
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
                'clim_dist_all_MMRR_results_%ikmrad.shp' % neigh_rad,
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

fig = plt.figure(figsize=(26,7))
gs = fig.add_gridspec(nrows=100, ncols=300)

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
emp_pval = np.mean(np.array(null_fx) >= fx_size)
emp_pval_str = '\n\nEMPIRICAL P-VALUE: %0.3f\n' % emp_pval
print(emp_pval_str)


#........................
# PART C: empirical p-val
#........................

# plot clim_dist effect size vs. empirical null dist
ax = fig.add_subplot(gs[80:, 215:])
ax.hist(null_fx, bins=100, alpha=0.7)
ax.axvline(fx_size, *ax.get_ylim(), color='red', linewidth=2)
ax.set_xlabel(('$\\Delta \\beta_{clim\\_dist}\\ /\\ '
               '\\Delta\\overline{x}_{lat}$'),
              fontdict=axlabel_fontdict)
ax.set_ylabel('freq', fontdict=axlabel_fontdict)
xticks = np.linspace(ax.get_xticks()[0], ax.get_xticks()[-1], 5) 
ax.set_xticks(xticks, ['%0.3f' % t for t in xticks])
ax.set_yticks(())
ax.set_ylabel('')
ax.tick_params(labelsize=ticklabel_size)
ax.text(0.6*ax.get_xlim()[1],
        0.8*ax.get_ylim()[1],
        '%0.3f' % fx_size,
        color='red',
        size=12,
       )
ax.text(0.6*ax.get_xlim()[1],
        0.6*ax.get_ylim()[1],
        '$(P=%0.3f)$' % emp_pval,
        color='red',
        size=12,
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
sns.regplot(x=np.abs(all_loop_MMRR_res_gdf['mean_lat']),
            y=all_loop_MMRR_res_gdf['clim_dist'],
            scatter=False,
            ax=ax,
           )
ax.set_xlabel('$|\\overline{x}_{lat}|$', fontdict=axlabel_fontdict)
ax.set_ylabel('$\\beta_{clim\\_dist}$', fontdict=axlabel_fontdict)
ax.tick_params(labelsize=ticklabel_size)

# add label for part B
ax.text(ax.get_xlim()[0] - (0.2*(ax.get_xlim()[1]-ax.get_xlim()[0])),
        1.05*ax.get_ylim()[1],
        'B.', size=partlabel_size, weight='bold')



#........................
# PART A: summary hex map
#........................

# NOTE: based on code taken from:
    # https://towardsdatascience.com/uber-h3-for-data-analysis-with-python-1e54acdcc908
ax = fig.add_subplot(gs[:, :190])

# add bottom axes for a colorbar
divider = make_axes_locatable(ax)
bcax = divider.append_axes('bottom', size='7%', pad=0.2)

# make dataframe to hold h3-converted data
h3_df = pd.DataFrame([], columns=['row_id', 'h3_id',
                                  'h3_geo_boundary', 'h3_centroid'])

# loop over results rows and convert to H3 hexes
for i, row in all_loop_MMRR_res_gdf.iterrows():
    p = row['geometry']
    if isinstance(p, MultiPolygon):
        ps = list(p)
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
for i, row in h3_gdf.iterrows():
    mean_clim_dist_val = float(all_loop_MMRR_res_gdf[
        all_loop_MMRR_res_gdf.intersects(row['geometry'])].mean()['clim_dist'])
    mean_clim_dist_vals.append(mean_clim_dist_val)
assert len(mean_clim_dist_vals) == len(h3_gdf)
h3_gdf['clim_dist_mean'] = mean_clim_dist_vals

# transform to equal-area projection and plot
subnational.plot(color='none',
                 edgecolor='black',
                 zorder=0,
                 ax=ax,
                 alpha=0.6,
                )
countries.plot(color='none',
                            edgecolor='black',
                            linewidth=1,
                            zorder=1,
                            ax=ax,
                            )
hex_subplot = h3_gdf.to_crs(8857).plot('clim_dist_mean',
                                       cmap='magma',
                                       alpha=0.9,
                                       zorder=2,
                                       ax=ax,
                                      )
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_xticks(())
ax.set_yticks(())
ax.set_title('')

# clip off longitudinally and latitudinally
ax.set_xlim(0.85 * ax.get_xlim()[0], 0.87*ax.get_xlim()[1])

ax.set_ylim(1.11*np.min(h3_gdf.to_crs(8857).geometry.centroid.y),
            1.175*np.max(h3_gdf.to_crs(8857).geometry.centroid.y))

# add label for part A
ax.text(1.08*ax.get_xlim()[0],
        1.065*ax.get_ylim()[1],
        'A.', size=partlabel_size, weight='bold')

# add the colorbar
# NOTE: need to create a custom ScalarMappable to feed into the colormap call
sm = plt.cm.ScalarMappable(cmap='magma',
                        norm=plt.Normalize(vmin=np.min(h3_gdf.clim_dist_mean),
                                           vmax=np.max(h3_gdf.clim_dist_mean)))
plt.colorbar(sm, cax=bcax, orientation='horizontal')
xticks = np.linspace(np.min(h3_gdf.clim_dist_mean),
                     np.max(h3_gdf.clim_dist_mean), 5)
bcax.set_xlabel('$average\\ \\beta_{clim\\_dist}$',
                fontdict=axlabel_fontdict)
bcax.set_xticks(xticks, ['%0.2f' % t for t in xticks], size=ticklabel_size)
bcax.set_ylabel('')
bcax.set_yticks(())



# adjust subplots and write to disk
fig.subplots_adjust(top=0.92,
                    bottom=0.12,
                    right=0.98,
                    left=0.03,
                    wspace=0.5,
                    hspace=0.3,
                   )

fig.savefig('FIG5_seas_dist_vs_clim_dist.png', dpi=700)


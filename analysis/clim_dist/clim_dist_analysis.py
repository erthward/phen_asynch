#!/usr/bin/env python
# coding: utf-8

################################################################################
# TODO:

    # need to standardize climate dists (and all other vals) within each MMRR
    # model? (otherwise, if clim dist is characteristically lower in tropics
    # than clim_dist coeff should also be smaller, but that reflects nothing
    # meaningful!)

################################################################################

# py packages
import geopandas as gpd
import numpy as np
import glob
import json
import time
import os
import sys
import re
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
import itertools

# local modules
import helper_fns as hf
from MMRR import MMRR


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
# plot Whittaker biomes?
plot_whittaker = True
# save CSVs of results?
save_all_results = False


#####################################
# PART I: CLUSTER HIGH-ASYNCH REGIONS
#####################################

# load the country boundaries (to use for random point generation)
cntry = gpd.read_file(hf.COUNTRIES_DATA_DIR + 'countries.shp')

# load asynch data
asynch = rxr.open_rasterio('../../results/maps/NIRv_global_asynch.tif')[2]

# determine if in interactive mode
interactive = hasattr(sys, 'ps1')
if interactive:
    loop_vals = itertools.product(dbscan_eps_vals,
                                  dbscan_minsamp_vals,
                                  alpha_vals)
else:
    loop_vals = [[float(arg.strip()) for arg in sys.argv[1:]]]

# loop analysis over param vals
for (dbscan_eps, dbscan_minsamp, alpha) in loop_vals:

    # wrap whole thing in a try statement, so that loop can skip over
    # parameterizations that produce less than min_n_clusts
    try:

        print(('#'*80+'\n')*4)
        print('dbscan_eps: %0.3f\n\n' % dbscan_eps)
        print('dbscan_minsamp: %0.3f\n\n' % dbscan_minsamp)
        print('alpha: %0.3f\n\n' % alpha)

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
                min_samples=dbscan_minsamp * maxes.shape[0]).fit(coords)
        # NOTE: ADD 1 TO CLUSTER LABELS, SO THAT FIRST CLUSTER == 1, NOT 0
        db.labels_ += 1
        # get clusters' mean latitudes (to be able to plot in order, N to S
        clust_mean_lats = {}
        for clust_i in np.unique(db.labels_):
            if clust_i > 0:
                mean_lat = np.mean(coords[db.labels_==clust_i,1])
                clust_mean_lats[clust_i] = mean_lat

        # check number of clusters identified
        print('\n\n%i CLUSTERS IDENTIFIED\n\n' % len(clust_mean_lats))
        assert len(clust_mean_lats) >= min_n_clusts

        # remap cluster labels so that they go in latitudinal order (N->S)
        NtoS_clust_labels = sorted(clust_mean_lats.keys(),
                                   key=lambda k: clust_mean_lats[k])[::-1]
        new_labels_dict = {0:0}
        for i, new_label in enumerate(NtoS_clust_labels):
            new_labels_dict[new_label] = i+1
        new_labels_list = []
        for label in db.labels_:
            new_labels_list.append(new_labels_dict[label])
        # make map
        fig = plt.figure(figsize=(22,5))
        gs = fig.add_gridspec(10, 3, width_ratios=[0.5,0.25,0.25])
        ax_map = fig.add_subplot(gs[:,0])
        cntry.to_crs(4326).plot(color='none', edgecolor='k', linewidth=0.25, ax=ax_map)
        ax_map.scatter(coords[:, 0], coords[:, 1], c=db.labels_, cmap='tab20', s=1)
        ax_map.scatter(coords[:, 0], coords[:, 1], c=db.labels_>0, cmap='Greys_r',
                       alpha=db.labels_==0, s=1)

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
                # coerce to MultiPolygon, if necessary
                if isinstance(alpha_hull, MultiPolygon):
                    pass
                else:
                    alpha_hull = MultiPolygon([alpha_hull])
                # save the MultiPolygon
                regs.append(alpha_hull)
                # add all constituent polygons to map
                for p in alpha_hull.geoms:
                    ax_map.add_patch(PolygonPatch(p, alpha=0.2, color='black'))
                cent = alpha_hull.centroid.coords.xy
                ax_map.text(cent[0][0], cent[1][0], str(clust_i), color='black',
                            size=14, weight='bold', alpha=0.85)
        # coerce all to MultiPolygons (to be able to use for raster clipping
        regs = [reg if isinstance(reg, MultiPolygon) else MultiPolygon([reg]) for reg in regs]


        ####################################
        # PART II: DRAW PTS AND RUN ANALYSIS
        ####################################

        regs_pts = [hf.generate_random_points_in_polygon(n_pts, reg) for reg in regs]

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
        geo_dist_colm = []
        reg_colm = []
        x1_colm = []
        x2_colm = []
        y1_colm = []
        y2_colm = []

        if run_MMRR:
            MMRR_res = {}

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
        mean_clim_dists = []
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

            # drop clim dists for points without ts dists, and vice versa
            not_missing = np.where(np.nansum(seas_dist, axis=0)>0)[0]
            seas_dist = seas_dist[:, not_missing][not_missing,:]
            clim_dist = clim_dist[:, not_missing][not_missing,:]
            pts = pts[not_missing, :]
            still_not_missing = np.where(np.nansum(clim_dist, axis=0)>0)[0]
            seas_dist = seas_dist[:, still_not_missing][still_not_missing,:]
            clim_dist = clim_dist[:, still_not_missing][still_not_missing,:]
            pts = pts[still_not_missing, :]

            print(('\n%i points remain after dropping locations without seasonality '
                  'data\n' % seas_dist.shape[0]))
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
                # run model
                res = MMRR(seas_dist, [hf.standardize_array(clim_dist), geo_dist],
                           ['clim_dist', 'geo_dist'])
                # print and store results
                for k, v in res.items():
                    print('\n\t%s: %s' % (k, str(v)))
                MMRR_res[reg_name] = res

                # extract the lower triangular values and scatter them
                indices = np.tril_indices(seas_dist.shape[0])

                seas_dist_vals = seas_dist[indices]
                clim_dist_vals = clim_dist[indices]
                geo_dist_vals = geo_dist[indices]

                # store the dists
                dist_dict[reg_name] = {'clim': clim_dist_vals,
                                       'seas': seas_dist_vals,
                                      }

                # add to DataFrame columns
                seas_dist_colm.extend(seas_dist_vals)
                clim_dist_colm.extend(clim_dist_vals)
                reg_colm.extend([reg_name]*len(seas_dist_vals))
                x1_colm.extend(pts[indices[0],0])
                x2_colm.extend(pts[indices[1],0])
                y1_colm.extend(pts[indices[0],1])
                y2_colm.extend(pts[indices[1],1])
                geo_dist_colm.extend(geo_dist_vals)

                # store the region's mean latitude and mean ppt
                mean_lat = [*reg_poly.centroid.coords.xy[1]][0]
                mean_ppt = float(np.mean(ppt.rio.clip(reg_poly)).values)
                mean_tmp = float(np.mean(tmp.rio.clip(reg_poly)).values)
                mean_clim_dist = np.mean(clim_dist_vals)
                mean_lats.append(mean_lat)
                mean_ppts.append(mean_ppt)
                mean_tmps.append(mean_tmp)
                mean_clim_dists.append(mean_clim_dist)

        # make results df

        print('\n\n\tPOINT 0\n\n')

        df = pd.DataFrame({'seas_dist': seas_dist_colm,
                           'clim_dist': clim_dist_colm,
                           'geo_dist': geo_dist_colm,
                           'reg': reg_colm,
                           'x1': x1_colm,
                           'x2': x2_colm,
                           'y1': y1_colm,
                           'y2': y2_colm,
                          })

        print('\n\n\tPOINT 1\n\n')

        # write results to disk
        if save_all_results:
            df.to_csv('clim_dist_results_multiregion%s_eps%0.2f_minsamp%0.2f_alpha*0.2f.csv'
                      % (file_suffix, dbscan_eps, dbscan_minsamp, alpha), index=False)

        print('\n\n\tPOINT 2\n\n')

        if run_MMRR:
            MMRR_res_df = pd.DataFrame(MMRR_res).T.reset_index()
            MMRR_res_df.columns = ['region']+[*MMRR_res_df.columns][1:]
            # add the mean latitude and mean ppt columns
            MMRR_res_df['mean_lat'] = mean_lats
            MMRR_res_df['mean_ppt'] = mean_ppts
            MMRR_res_df['mean_tmp'] = mean_tmps
            MMRR_res_df['mean_clim_dist'] = mean_tmps
            if save_all_results:
                MMRR_res_df.to_csv('clim_dist_MMRR_results_multiregion%s_eps%0.2f_minsamp%0.2f_alpha*0.2f.csv'
                                   % (file_suffix, dbscan_eps, dbscan_minsamp, alpha), index=False)


            # plot clim_dist MMRR coeff as a function of mean lat and mean ppt
            ax_scat = fig.add_subplot(gs[:5, 1])

            ax_scat.scatter(np.abs(MMRR_res_df['mean_lat']),
                            MMRR_res_df['clim_dist'],
                            c='black',
                            s=40,
                            alpha=0.9,
                            )

            sns.regplot(x=np.abs(MMRR_res_df['mean_lat']),
                        y=MMRR_res_df['clim_dist'],
                        scatter=False,
                        ax=ax_scat,
                       )

            for r,x,y in zip(MMRR_res_df['region'].values,
                             np.abs(MMRR_res_df['mean_lat'].values),
                             MMRR_res_df['clim_dist'].values):
                ax_scat.text(x, y, r, size=16, alpha=0.75, color=reg_cols[int(r)-1])

            ax_scat.set_xlabel('absolute mean cluster latitude (|degrees|)',
                               fontdict={'fontsize': 14})
            ax_scat.set_ylabel('MMRR standardized climate\ndistance coefficient',
                               fontdict={'fontsize': 14})
            ax_scat.tick_params(labelsize=11)


            # plot mean clim_dist against mean lat (as a 'null')
            ax_scat_null = fig.add_subplot(gs[5:, 1])
            ax_scat_null.scatter(np.abs(MMRR_res_df['mean_lat']),
                            MMRR_res_df['mean_clim_dist'],
                            c='black',
                            s=40,
                            alpha=0.9,
                            )
            sns.regplot(x=np.abs(MMRR_res_df['mean_lat']),
                        y=MMRR_res_df['mean_clim_dist'],
                        scatter=False,
                        ax=ax_scat_null,
                       )
            for r,x,y in zip(MMRR_res_df['region'].values,
                             np.abs(MMRR_res_df['mean_lat'].values),
                             MMRR_res_df['mean_clim_dist'].values):
                ax_scat_null.text(x, y, r, size=16, alpha=0.75, color=reg_cols[int(r)-1])
            ax_scat_null.set_xlabel('absolute mean cluster latitude (|degrees|)',
                               fontdict={'fontsize': 14})
            ax_scat_null.set_ylabel('mean pairwise climate\ndistance of cluster',
                               fontdict={'fontsize': 14})
            ax_scat_null.tick_params(labelsize=11)



            # plot on biomes
            if plot_whittaker:
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
                ax_whit = fig.add_subplot(gs[:, 2])
                divider = make_axes_locatable(ax_whit)
                cax_whit = divider.append_axes('right', size='5%', pad=0.1)
                ax_whit.add_collection(p)

                for b,c in zip(biomes, centroids):
                    ax_whit.text(c[0], c[1], b)

                for r,x,y in zip(MMRR_res_df['region'].values,
                                 MMRR_res_df['mean_tmp'].values,
                                 MMRR_res_df['mean_ppt'].values):
                    ax_whit.text(x, y, r, size=16, alpha=0.5)

                scat = ax_whit.scatter(MMRR_res_df['mean_tmp'], MMRR_res_df['mean_ppt'],
                                   c=MMRR_res_df['clim_dist'],
                                   cmap='plasma_r', s=45, alpha=0.85, edgecolor='k',
                                   linewidth=1)

                cbar = plt.colorbar(scat, cax=cax_whit)
                cbar.set_label('MMRR climate-distance coefficient',
                               fontdict={'fontsize': 18})

                ax_whit.set_xlabel('MAT ($^{\circ}C$)',
                              fontdict={'fontsize': 18})
                ax_whit.set_ylabel('MAP ($mm$)',
                              fontdict={'fontsize': 18})
                ax_whit.tick_params(labelsize=14)

        fig.subplots_adjust(bottom=0.06,
                            top=0.94,
                            left=0.06,
                            right=0.94,
                            wspace=0.25,
                           )

        fig.savefig(os.path.join('asynch_clust_results',
                                 ('asynch_clusts_and_MMRR_results'
                                  '%s_eps%0.2f_minsamp%0.2f_'
                                  'alpha%0.2f.png') % (file_suffix,
                                                       dbscan_eps,
                                                       dbscan_minsamp,
                                                       alpha)), dpi=800)

    except Exception as e:
        print('\n\n\nException raised: %s\n\nMoving on...\n\n\n' % e)

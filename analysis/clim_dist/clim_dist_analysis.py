#!/usr/bin/env python
# coding: utf-8

# py packages
import geopandas as gpd
import numpy as np
import glob
import json
import time
import os
#from osgeo import gdal
import random
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
import seaborn as sns
import math
import itertools
import rasterio as rio
import matplotlib.pyplot as plt
from shapely.ops import unary_union
from shapely.geometry import Polygon, Point, Multipolygon
from scipy.spatial import ConvexHull
from collections import Counter as C

# local modules
import helper_fns as hf


# TODO:
    # perhaps I should just be trying to calculate neighborhood multivariate
    # variance in standardized bioclim variables as a RF covariate (rather
    # than this inter-point comparison thing)?

    # change this script to sample points in a stratified random way,
    # then assign regions in some data-driven way

    # account for/avoid spatial autocorrelation in some way?


# load the country boundaries (to use for random point generation)
cntry = gpd.read_file(hf.COUNTRIES_DATA_DIR + 'countries.shp')

n_pts = 1000  # Enter in the number greater than random points you need


# define regions for:
# high asynchrony:
    # temperate: CA and SW
ha_tmp = Polygon([[-123.588671875,42.1726562502269],
                  [-120.951953125,35.75664368702668],
                  [-109.43828125,31.96275411512188],
                  [-105.746875,32.03728956520954],
                  [-106.362109375,43.653059057103164],
                  [-123.588671875,42.1726562502269]
                 ])
    # tropical: N Andes and N SA
ha_trp = Polygon([[-76.84073042457382,-0.40960685149885684],
                  [-58.12002729957382,-0.40960685149885684],
                  [-58.12002729957382,11.466203915130043],
                  [-76.84073042457382,11.466203915130043],
                  [-76.84073042457382,-0.40960685149885684]
                 ])
#ha_trp = Polygon([[-76,8.5],
#                  [-76,-1],
#                  [-58,-1],
#                  [-58,8.5],
#                  [-76,8.5]
#                 ])
# low asynchrony:
    # temperate: US South and SE
la_tmp = Polygon([[-95.39198282257864,30.76200580067661],
                  [-81.68104532257864,30.76200580067661],
                  [-81.68104532257864,37.797946552183234],
                  [-95.39198282257864,37.797946552183234],
                  [-95.39198282257864,30.76200580067661]
                 ])
    # tropical: W Amazon
la_trp = Polygon([[-72.70987104957382,-11.49384470495583],
                  [-58.12002729957382,-11.49384470495583],
                  [-58.12002729957382,-1.9912406261155682],
                  [-72.70987104957382,-1.9912406261155682],
                  [-72.70987104957382,-11.49384470495583]
                 ])
regs = [ha_tmp, ha_trp, la_tmp, la_trp]

# draw random points for each region
regs_pts = [hf.generate_random_points_in_polygon(n_pts, reg) for reg in regs]

# create the figure
fig, ax = plt.subplots(1,1)
fig.suptitle(('seasonal distance vs. climatic distance, '
              'across latitude and asynchrony'), size=20)

# list of region names
reg_names = ['high asynch: temperate',
             'high asynch: tropical',
             'low asynch: temperate',
             'low asynch: tropical'
            ]

# create colors for the 4 regions
reg_cols = ["#7428ed", # ha_tmp
            "#18f04a", # ha_trp
            "#a6a4b3", # la_tmp
            "#8db08f", # la_trp
           ]

# dict to store the distances
dist_dict = {}

# get the nodata val
nodata_val = rio.open(hf.BIOCLIM_INFILEPATHS[0]).nodata

# columns for pandas DataFrame
seas_dist_colm = []
clim_dist_colm = []
reg_colm = []
x1_colm = []
x2_colm = []
y1_colm = []
y2_colm = []

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
    seas_dist = hf.get_raster_info_points(hf.COEFFS_FILE, pts, 'ts_pdm')

    # drop clim dists for points without ts dists, and vice versa
    not_missing = np.where(np.nansum(seas_dist, axis=0)>0)[0]
    seas_dist = seas_dist[:, not_missing][not_missing,:]
    clim_dist = clim_dist[:, not_missing][not_missing,:]
    still_not_missing = np.where(np.nansum(clim_dist, axis=0)>0)[0]
    seas_dist = seas_dist[:, still_not_missing][still_not_missing,:]
    clim_dist = clim_dist[:, still_not_missing][still_not_missing,:]

    print(('\n%i points remain after dropping locations without seasonality '
          'data\n' % seas_dist.shape[0]))

    # extract the lower triangular values and scatter them
    indices = np.tril_indices(seas_dist.shape[0])

    seas_dist_vals = seas_dist[indices]
    clim_dist_vals = clim_dist[indices]

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
                           'seas': seas_dist_vals
                          }

    # add to DataFrame columns
    seas_dist_colm.extend(seas_dist_vals)
    clim_dist_colm.extend(clim_dist_vals)
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
df.to_csv('clim_dist_results.csv', index=False)




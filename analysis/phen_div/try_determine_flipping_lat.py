#--------
# imports
#--------

from sklearn.linear_model import LinearRegression
from rasterio.fill import fillnodata
from geopy.distance import geodesic
from greatcircle import GreatCircle
from shapely.geometry import Point
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from copy import deepcopy
from pprint import pprint
import rasterstats as rs
import geopandas as gpd
import rioxarray as rxr
import numpy as np
import glob
import json
import time
import os

import matplotlib.pyplot as plt
from shapely.ops import unary_union
from shapely.geometry import Polygon, Point
from scipy.spatial import ConvexHull
from collections import Counter as C

DATA_DIR = '/home/deth/Desktop/CAL/research/projects/seasonality/results/maps/'
COEFFS_FILE = DATA_DIR + 'NIRv_global_coeffs.tif'

coeffs = rxr.open_rasterio(COEFFS_FILE)


#-----------
# set params
#-----------

# default missing-data val
DEFAULT_VAL = -9999.0

# names of the bands saved into the TFRecord files
INBANDS = ['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']



#-----------
# helper fns
#-----------


def make_design_matrix():
    """
    Makes and returns the regression's design matrix, a 365 x 5 numpy array
    in which the columns contain, in order:
        - 1s (for the constant);
        - sin and cos of annual-harmonic days of the year
          (i.e. days are expressed in radians from 0 to 2pi);
        - sin and cos of the semiannual-harmonic days of the year
          (i.e. days are expressed in radians from 0 to 4pi).
    """
    # get 1 year of daily values, expressed in radians, 1 rotation/yr
    annual_radian_days = np.linspace(0, 2*np.pi, 366)[:365]
    # get 1 year of daily values, expressed in radians, 2 rotations/yr
    semiannual_radian_days = np.linspace(0, 4*np.pi, 366)[:365] % (2 * np.pi)
    # get the harmonic values of those
    sin1 = np.sin(annual_radian_days)
    cos1 = np.cos(annual_radian_days)
    sin2 = np.sin(semiannual_radian_days)
    cos2 = np.cos(semiannual_radian_days)
    # add a vector of 1s for the constant term, then recast as a 365 x 5 array,
    # to use as the covariate values in the regression
    design_mat = np.array([np.ones(sin1.shape), sin1, cos1, sin2, cos2]).T
    return design_mat


def get_time_series(coeffs, design_mat, normalize=True):
    """
    Calculates the time series at the pixel the coeffs came from,
    using the coefficients for the
    constant and the sin and cosine terms of the annual and semiannual
    harmonic components. Returns the time series as a numpy array.
    """
    # multiply the pixel's set of coefficients by the design mat, then sum
    # all the regression terms
    # NOTE: the coeffs are a numpy array of shape (5,);
    #       the design matrix is a numpy array of shape (365, 5);
    #       the multiplication operates row by row, and elementwise on each
    #       row, returning a new (365, 5) array that, when summed across the
    #       columns (i.e. the regression's linear terms) gives a (365,) vector
    #       of fitted daily values for the pixel
    ts = np.sum(coeffs * design_mat, axis=1)
    if normalize:
        ts = (ts-np.min(ts))/(np.max(ts)-np.min(ts))
    return ts


#-----
# main
#-----

dm = make_design_matrix()

min_lat = np.min(coeffs.y.values)
max_lat = np.max(coeffs.y.values)

# TODO: REPLACE WITH LOOP OVER ALL LATS
for lon in [-65]:
    fig, ax = plt.subplots(1,1)
    coeffs_lon = coeffs.sel(x=lon, method='nearest')
    for lat_i, lat in enumerate(coeffs_lon.y.values):
        color_frac = (lat-min_lat)/(max_lat-min_lat)
        color = plt.cm.viridis(color_frac)
        coeffs_lon_lat = coeffs_lon[:, lat_i].values
        if np.nan not in coeffs_lon_lat:
            ts = get_time_series(coeffs_lon_lat, dm)
            ax.plot(ts, color=color, alpha=0.5)

    fig.show()



# NOTES:
    # IDEA:
    # for each lon:
        # start assuming 0 is the 'flip' point
        # for pairs of ts stepping away from 0 lat,
        # calculate improvement in R2 of real N ~ rotated S ts (vs both real);
        # develop metric of goodness of fit of flip point;
        # then use binary search moving the flip point away from 0 to optimize
        # that metric

    # once I've done that for lons across each continent then fit a spline
    # between to determine flipping-latitude for each longitude

    # then rotate time series at those points

    # then fit EOF


abs_lats = sorted(list(set(np.round(np.abs(coeffs.y.values), 3))))
step = abs_lats[1]-abs_lats[0]
flip_points = [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7,
               8, -8, 9, -9, 10, -10, 11, -11, 12, -12]
for lon in [-65]:
    coeffs_lon = coeffs.sel(x=lon, method='nearest')
    all_R2_diffs = {}
    all_mean_R2_diffs = {}
    for flip_point in flip_points:
        R2_diffs = []
        for i in range(700):
            lat_N = flip_point+((i+1)*step)
            lat_S = flip_point-((i+1)*step)
            coeffs_N = coeffs_lon.sel(y=lat_N, method='nearest').values
            coeffs_S = coeffs_lon.sel(y=lat_S, method='nearest').values
            if (np.sum(np.isnan(coeffs_N)) == 0 and
                np.sum(np.isnan(coeffs_S)) == 0):
                ts_N = get_time_series(coeffs_N, dm)
                ts_S = get_time_series(coeffs_S, dm)
                R2_raw = stats.pearsonr(ts_N, ts_S)[0]**2
                R2_rot = stats.pearsonr(ts_N, np.concatenate((ts_S[183:],
                                                         ts_S[:183])))[0]**2
                R2_diff = R2_rot-R2_raw
                R2_diffs.append(R2_diff)
        all_R2_diffs[flip_point] = R2_diffs
        all_mean_R2_diffs[flip_point] = np.mean(R2_diffs)
    print('Lat with max mean difference in R2s:')
    print([k for k, v in all_mean_R2_diffs.items() if v ==
           np.max([*all_mean_R2_diffs.values()])][0])



min_fp = np.min(flip_points)
max_fp = np.max(flip_points)
fig, ax = plt.subplots(1,1)
for fp, R2s in all_R2_diffs.items():
    color_frac = (fp-min_fp)/(max_fp-min_fp)
    color = plt.cm.viridis(color_frac)
    ax.plot(R2s, color=color, alpha=0.5)
fig.show()

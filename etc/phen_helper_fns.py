#!/bin/python
# phen_helper_fns.py

"""
Functions to help with downstream analyses
of the asynch and harmonic-regression maps.
"""

#--------
# imports
#--------

from sklearn.linear_model import LinearRegression
from rasterio.fill import fillnodata
from geopy.distance import geodesic
from shapely.geometry import Point
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from copy import deepcopy
from pprint import pprint
import rasterstats as rs
import rioxarray as rxr
import geopandas as gpd
import rasterio as rio
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

# set the neighborhood radii (in km)
neigh_rads = [50, 100, 150]

# get the variables' filepaths
DATA_DIR = ('/home/deth/Desktop/CAL/research/projects/seasonality/'
            'seasonal_asynchrony/data')
EXTERNAL_DATA_DIR = '/media/deth/SLAB/diss/3-phn/GEE_outputs/final/'
EXTERNAL_RF_DATA_DIR = '/media/deth/SLAB/diss/3-phn/rf_data'
COEFFS_FILE = os.path.join(EXTERNAL_DATA_DIR, 'NIRv_coeffs.tif')
COEFFS_STRICT_FILE = os.path.join(EXTERNAL_DATA_DIR, 'NIRv_STRICT_coeffs.tif')
ASYNCH_FILES = {rad: os.path.join(EXTERNAL_DATA_DIR,
                'NIRv_STRICT_asynch_%ikm.tif' % rad) for rad in neigh_rads}
BOUNDS_DIR = os.path.join(DATA_DIR, 'bounds')
BIOCLIM_DIR = os.path.join(DATA_DIR, 'bioclim')
BIOCLIM_INFILEPATHS = glob.glob(os.path.join(BIOCLIM_DIR,"wc2.1_2.5m_bio_*.tif"))
BIOCLIM_INFILEPATHS = sorted(BIOCLIM_INFILEPATHS)


#-----------
# set params
#-----------

# pattern that occurs just before the file number in each input file's name
PATT_B4_FILENUM = '-OUT-'

# kernel size used by GEE to output the TFRecord files
KERNEL_SIZE = 60

# default missing-data val
DEFAULT_VAL = -9999.0

# names of the bands saved into the TFRecord files
INBANDS = ['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']
OUTBANDS = ['asynch', 'asynch_R2', 'asynch_euc', 'asynch_euc_R2', 'asynch_n']


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


def calc_time_series(design_mat, coeffs=None, i=None, j=None, patch=None):
    """
    Calculates the time series at pixel i,j, using the coefficients for the
    constant and the sin and cosine terms of the annual and semiannual
    harmonic components. Returns the time series as a numpy array.
    """
    assert ((i is not None and j is not None and patch is not None and coeffs is None) or
            (coeffs is not None and i is None and j is None and patch is None))
    # multiply the pixel's set of coefficients by the design mat, then sum
    # all the regression terms
    # NOTE: the coeffs are a numpy array of shape (5,);
    #       the design matrix is a numpy array of shape (365, 5);
    #       the multiplication operates row by row, and elementwise on each
    #       row, returning a new (365, 5) array that, when summed across the
    #       columns (i.e. the regression's linear terms) gives a (365,) vector
    #       of fitted daily values for the pixel
    if patch is None:
        ts = np.sum(coeffs * design_mat, axis=1)
    else:
        ts = np.sum(patch[:, i, j] * design_mat, axis=1)
    return ts


def calc_euc_dist(a1, a2):
    """
    Calculates the Euclidean distance between two 1d, length-n numpy arrays.

    Returns the distance as a float.
    """
    dist = np.sqrt(np.sum((a1 - a2)**2))
    return dist


def standardize_array(a):
    """
    Returns a standardized version of the input array
    """
    return (a - np.mean(a))/np.std(a)


def minmax_rescale_array(a):
    """
    Returns a min-max rescaled (i.e., min-max normalized) version of the input array
    """
    return (a-np.min(a))/(np.max(a)-np.min(a))


def get_tolerancefilled_arr(arr, tol):
    """
    uses rasterio.fill.fillnodata to interpolate missing values
    in an array, but only for nan-patches that are < 2*tol cells wide
    in both the horizontal and vertical directions
    """
    # get binary masking array, where 0s are nans, 1s are not
    mask = np.int8(np.invert(np.isnan(arr)))
    # loop over cells and sum all cells within tol linear dist in all 4 dirs,
    # then use sum to determine if a nan cell has non-nans within that nhood
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            # only compute this for missing cells
            if mask[i,j] == 0:
                # get the window
                wind_x_min = max(j-tol, 0)
                wind_x_max = min(j+tol, arr.shape[1])
                wind_y_min = max(i-tol, 0)
                wind_y_max = min(i+tol, arr.shape[0])
                # get count of non-nans within window
                tot = np.sum(np.invert(np.isnan(arr[wind_y_min:wind_y_max,
                                                    wind_x_min:wind_x_max])))
                # update mask so that this cell isn't filled,
                # if too far from non-nan
                if tot == 0:
                    mask[i,j] = -1
    # create a deepcopy (because fillnodata appears to update its arr in-place,
    # even though the docs seem to imply to the contrary)
    out_arr = deepcopy(arr)
    # fill and return the copy
    filled = fillnodata(out_arr, mask)
    return filled


def gapfill_and_rewrite_raster(rast_filepath, fill_tol=5):
    # read in the raster file
    f = rio.open(rast_filepath)
    rast = f.read()
    #rast = rxr.open_rasterio(rast_filepath, masked=True)

    # fill missing values in each band of the raster,
    # using rasterio.fill.fillnodata
    for n in range(rast.shape[0]):
        rast[n, :, :] = get_tolerancefilled_arr(rast[n, :, :],
                                                tol=fill_tol)

    # write new file
    out_profile = f.profile
    path_and_basename, ext = os.path.splitext(rast_filepath)
    with rio.open(path_and_basename + '_FILLED_tol%i' % fill_tol + ext,
                       'w',
                       **f.profile) as dst:
        dst.write(rast.astype(out_profile['dtype']),
                  range(1, f.profile['count'] + 1))


def get_raster_info_points(rast_filepath, pts, return_format='vals',
                           standardize=False, fill_nans=False, fill_tol=5):
    """
    takes a raster (e.g. coeffs from the harmonic seasonality regression)
    and an nx2 np.array of n points' x and y (i.e. lon, lat) coordinates,
    returns either the values at those points (if return_format == 'vals';
    default), the fitted time series at those
    points (if return_format == 'ts'; ONLY VALID FOR COEFF RASTERS!),
    the pairwise distance matrix (if return_format
    == 'pdm'), of the pairwise seasonal distance matrix (if
    return_format=='ts_pdm'; ONLY VALID FOR COEFF RASTERS!)
    """
    # read in the raster file
    f = rio.open(rast_filepath)
    rast = f.read()
    #rast = rxr.open_rasterio(rast_filepath, masked=True)

    # fill missing values in each band of the raster,
    # using rasterio.fill.fillnodata, if requested
    if fill_nans:
        for n in range(rast.shape[0]):
            rast[n, :, :] = get_tolerancefilled_arr(rast[n, :, :],
                                                    tol=fill_tol)

    # get the cells' max lon- and lat-bound values
    #cell_max_lons, cell_max_lats = get_cell_lonlat_bounds(f)

    if 'ts' in return_format:
        # get the regression's design matrix
        design_mat = make_design_matrix()

        # matrix-multiply design_mat x coeffs to get
        # fitted time series for all pts
        ts_mat = np.zeros((pts.shape[0], 365)) * np.nan

    else:
        vals_mat = np.zeros((pts.shape[0], f.count)) * np.nan


    for row_i in range(pts.shape[0]):
        # get the array coords of the cell the point falls in
        #pt_cell_i, pt_cell_j = get_cell_coords_for_pt(lon=pts[row_i,0],
        #                                              lat=pts[row_i,1],
        #                                            cell_max_lons=cell_max_lons,
        #                                            cell_max_lats=cell_max_lats)
        pt_cell_i, pt_cell_j = f.index(pts[row_i,0], pts[row_i,1])

        # calculate the ts, if needed
        if 'ts' in return_format:
            ts = np.sum(rast[:, pt_cell_i, pt_cell_j] * design_mat,
                        axis=1)
            # standardize the ts?
            if standardize:
                ts = standardize_array(ts)
            # store the ts
            ts_mat[row_i, :] = ts

        # otherwise, just extract the raster's values
        else:
            vals = rast[:, pt_cell_i, pt_cell_j]
            vals_mat[row_i, :] = vals

    if 'pdm' in return_format:
        # calculate pairwise Euc dist matrix
        pw_dist_mat = np.zeros((pts.shape[0], pts.shape[0])) * np.nan
        for row_i in range(pw_dist_mat.shape[0]):
            for col_j in range(pw_dist_mat.shape[1]):
                if row_i == col_j:
                    dist = 0
                else:
                    if return_format == 'pdm':
                        dist = calc_euc_dist(vals_mat[row_i,:],
                                             vals_mat[col_j,:])

                    elif return_format == 'ts_pdm':
                        dist = calc_euc_dist(ts_mat[row_i,:], ts_mat[col_j, :])
                pw_dist_mat[row_i, col_j] = dist
        return pw_dist_mat
    else:
        # return the time series matrix, if pairwise distances not requested
        if return_format == 'ts':
            return ts_mat
        else:
            return vals_mat


def get_cell_lonlat_bounds(rio_file):
    """
    takes a rasterio opened-file object,
    returns a tuple of two arrays, the first for lons, the second for lats,
    each of which contains a cell's max coordinate values for each col (lons)
    or row (lats) in the raster
    """
    # get the res, min, and length values for lon and lat dimensions
    aff = rio_file.transform
    x_res = aff.a
    x_min = aff.c
    x_len = rio_file.width
    y_res = aff.e
    y_min = aff.f
    y_len = rio_file.height
    # get arrays of the max (i.e. east and south) lon and lat cell boundaries
    # for all the cells in the raster
    lons_east_cell_bounds = np.linspace(start=x_min+x_res,
                                        stop=(x_res*x_len)+(x_min+x_res),
                                        num=x_len)
    lats_south_cell_bounds = np.linspace(start=y_min+y_res,
                                         stop=(y_res*y_len)+(y_min+y_res),
                                         num=y_len)
    return lons_east_cell_bounds, lats_south_cell_bounds


def get_cell_coords_for_pt(lon, lat, cell_max_lons, cell_max_lats):
    """
    uses provided dictionary to crosswalk a point's coordinates
    with its cell's row,col (i.e. i,j) coordinates in the numpy array
    holding a raster's data
    """
    #londiff_gtz = (cell_max_lons - lon) > 0
    #latdiff_gtz = (cell_max_lats - lat) > 0
    #j = np.argmax(~londiff_gtz[:-1] == londiff_gtz[1:]) + 1
    #i = np.argmax(~latdiff_gtz[:-1] == latdiff_gtz[1:]) + 1
    j = np.where((cell_max_lons - lon) > 0)[0][0]
    i = np.where((cell_max_lats - lat) < 0)[0][0]
    return i,j


# generate n random points within a given polygon
def generate_random_points_in_polygon(n, polygon):
    points = []
    minx, miny, maxx, maxy = polygon.bounds
    while len(points) < n:
        pnt = Point(np.random.uniform(minx, maxx),
                    np.random.uniform(miny, maxy))
        if polygon.contains(pnt):
            points.append(pnt)
    return points


def calc_pw_clim_dist_mat(pts, nodata_val=-3.4e+38):
    """
    Calculates the pw bioclim-dist matrix for an nx2 numpy array
    containing the (lon, lat) coordinates for n points.
    """
    # TODO:
        # IMPLEMENT ARGUMENT TO NORMALIZE BIOCLIM DATA?

    # list to hold arrays of all bioclim vars' vals
    bioclim_vals_arrs = []
    for fn in BIOCLIM_INFILEPATHS:
        # extract vals for each file
        bioclim_vals = get_raster_info_points(fn, pts, 'vals')
        # mask out missing vals
        bioclim_vals[np.isclose(bioclim_vals, nodata_val)] = np.nan
        # standardize the vals
        stand_vals = (bioclim_vals -
                      np.nanmean(bioclim_vals))/np.nanstd(bioclim_vals)
        bioclim_vals_arrs.append(stand_vals)
    # concatenate into one array
    bioclim = np.concatenate(bioclim_vals_arrs, axis=1)
    # calculate pairwise dists
    clim_dist = np.zeros([pts.shape[0]]*2) * np.nan
    for i in range(pts.shape[0]):
        for j in range(i, pts.shape[0]):
            if i == j:
                clim_dist[i,j] = 0
            else:
                dist = calc_euc_dist(bioclim[i,:],
                                        bioclim[j,:])
                clim_dist[i,j] = dist
                clim_dist[j,i] = dist
    return clim_dist



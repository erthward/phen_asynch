#!/bin/python
# helper_fns.py

"""
Functions to help with downstream analyses
of the asynch and harmonic-regression maps.
"""

#--------
# imports
#--------

from sklearn.linear_model import LinearRegression
from geopy.distance import geodesic
from shapely.geometry import Point
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from pprint import pprint
import rasterstats as rs
import geopandas as gpd
#import tensorflow as tf
import rasterio as rio
import numpy as np
import glob
import json
import time
import os


#-----------
# set params
#-----------

# directory where the data and mixerfile live
if os.path.abspath('.').split('/')[1] == 'home':
    DATA_DIR = ('/home/deth/Desktop/stuff/berk/research/projects/seasonality/'
                'GEE_output')
else:
    DATA_DIR = '/global/home/users/drewhart/seasonality/GEE_output/SIF/'

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


def calc_time_series(patch, i, j, design_mat):
    """
    Calculates the time series at pixel i,j, using the coefficients for the
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
    ts = np.sum(patch[:, i, j] * design_mat, axis=1)
    return ts


def calc_euc_dist(a1, a2):
    """
    Calculates the Euclidean distance between two 1d, length-n numpy arrays.

    Returns the distance as a float.
    """
    dist = np.sqrt(np.sum((a1 - a2)**2))
    return dist


def get_raster_info_points(rast_filepath, pts, return_format='vals'):
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



import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from copy import deepcopy
import numpy as np
import scipy.stats
import xarray as xr
import rasterio as rio
import rioxarray as rxr
from datetime import datetime
import re, os, sys, time

delete_after_finished = True

# choose masking mode
masking_mode = 'default'
masking_suffix = '_STRICT' * (masking_mode == 'strict')

# data directory
coeffs_data_dir = '/global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives/'


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


def calc_r2(a1, a2):
    """
    Calculates the R^2 value between both input arrays
    """
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(a1, a2)
    return r_value**2


def calc_euc_dist(a1, a2):
    """
    Calculates the Euclidean distance between two 1d, length-n numpy arrays.

    Returns the distance as a float.
    """
    dist = np.sqrt(np.sum((a1 - a2)**2))
    return dist


def rescale_data(data, lo=0, hi=1):
    assert hi>lo, 'hi must be > lo!'
    norm_data = (data - np.min(data))/(np.max(data)-np.min(data))
    norm_data = lo + (norm_data*(hi-lo))
    return norm_data


def load_rs_coeffs(rs_coeffs_tif):
    """
    load the RS-based seasonality coeffs file, to be validated
    """
    rs_coeffs = rxr.open_rasterio(rs_coeffs_tif, masked=True)
    # NOTE: already in WGS84 latlon, and I'll just be extracting using
    #       latlon, so no need to transform
    return rs_coeffs


def get_fitted_seas(coeffs_rast, x, y, design_mat):
    """
    Calculates the predicted time series at pixel at x,y in a rioxarray raster,
    using the coefficients for the constant and the
    sin and cosine terms of the annual and semiannual
    harmonic components. Returns the time series as a numpy array.
    """
    # get the nearest pixel's coefficients from the rioxarray dataset
    coeffs = coeffs_rast[:, i, j].values

    # multiply the pixel's set of coefficients by the design mat, then sum
    # all the regression terms
    # NOTE: the coeffs are a numpy array of shape (5,);
    #       the design matrix is a numpy array of shape (365, 5);
    #       the multiplication operates row by row, and elementwise on each
    #       row, returning a new (365, 5) array that, when summed across the
    #       columns (i.e. the regression's linear terms) gives a (365,) vector
    #       of fitted daily values for the pixel
    pred = np.sum(coeffs * design_mat, axis=1)

    # min-max rescale it
    pred = rescale_data(pred)

    return pred


# load the RS-based coefficients files
rs_coeffs_tif_filename = os.path.join(coeffs_data_dir, '%s/%s_coeffs%s.tif')
nirv_filename = rs_coeffs_tif_filename % ('NIRv', 'NIRv', masking_suffix)
sif_filename = rs_coeffs_tif_filename % ('SIF', 'SIF', masking_suffix)
nirv = load_rs_coeffs(nirv_filename)
sif = load_rs_coeffs(sif_filename)
# light check on coregistration, since I already know by now that they
# basically need to be coregistered
assert nirv.rio.crs == sif.rio.crs and np.all(nirv.shape == sif.shape)

# make design matrix used to estimate rs-based coeffs' fitted seasonality
design_mat = make_design_matrix()

# dataset to store results
r2s = deepcopy(nirv[0]*0)

# loop over pixels, calclating and storing R^2 between
# fitted NIRv and SIF seasonal phenologies
for i in range(r2s.shape[0]):
    print(f'\n\tprocessing row {i}...\n')
    for j in range(r2s.shape[1]):

        if np.nan in nirv[:,i,j] or np.nan in sif[:,i,j]:
            r2s[i,j] = np.nan

        else:
            nirv_preds = get_fitted_seas(nirv, i, j, design_mat)
            sif_preds = get_fitted_seas(sif, i, j, design_mat)
            assert len(nirv_preds) == len(sif_preds) == 365
            r2 = calc_r2(nirv_preds, sif_preds)
            r2s[i, j] = r2

# write results to disk
r2s.rio.to_raster(os.path.join(coeffs_data_dir, 'NIRv_SIF_phen_R2s.tif'))


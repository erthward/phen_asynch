#--------
# imports
#--------

import numpy as np
import rasterio as rio
import matplotlib.pyplot as plt
from eofs.standard import Eof
import cartopy.crs as ccrs

from sklearn.linear_model import LinearRegression
from geopy.distance import geodesic
from shapely.geometry import Point
from scipy.spatial import cKDTree
from pprint import pprint
import rasterstats as rs
import geopandas as gpd
import glob
import json
import time
import os



def get_lats_lons(transform, dims):
    """
    get lat and lon values for a raster with the given geotransform and dims
    """
    lat_dim = dims[0]
    lon_dim = dims[1]

    lat_min = transform[3]
    lon_min = transform[0]

    lat_res = transform[5]
    lon_res = transform[1]

    # NOTE: add half res to each min value, to indicate cell centers (rather
    # than the cell corners they represent by default)
    lat_min = lat_min + 0.5*lat_res
    lon_min = lon_min + 0.5*lon_res

    lat_max = lat_min + lat_dim*lat_res
    lon_max = lon_min + lon_dim*lon_res

    lats = np.linspace(lat_min, lat_max, lat_dim)
    lons = np.linspace(lon_min, lon_max, lon_dim)

    return lats, lons


def make_design_matrix(thin=False, thin_frac=0.1):
    """
    Makes and returns the regression's design matrix, a 365 x 5 numpy array
    in which the columns contain, in order:
        - 1s (for the constant);
        - sin and cos of annual-harmonic days of the year
          (i.e. days are expressed in radians from 0 to 2pi);
        - sin and cos of the semiannual-harmonic days of the year
          (i.e. days are expressed in radians from 0 to 4pi).
    """
    n_steps = 365
    if thin:
        n_steps = int(n_steps * thin_frac)
    # get 1 year of values, expressed in radians, 1 rotation/yr
    # NOTE: will be daily, unless thinned
    annual_radian_days = np.linspace(0, 2*np.pi, n_steps + 1)[:n_steps]
    # get 1 year of values, expressed in radians, 2 rotations/yr
    # NOTE: will be daily, unless thinned
    semiannual_radian_days = np.linspace(0, 4*np.pi,
                                         n_steps + 1)[:n_steps] % (2 * np.pi)
    # get the harmonic values of those
    sin1 = np.sin(annual_radian_days)
    cos1 = np.cos(annual_radian_days)
    sin2 = np.sin(semiannual_radian_days)
    cos2 = np.cos(semiannual_radian_days)
    # add a vector of 1s for the constant term, then recast as a 365 x 5 array,
    # to use as the covariate values in the regression
    design_mat = np.array([np.ones(sin1.shape), sin1, cos1, sin2, cos2]).T
    return design_mat


def calc_time_series_cube(rast, design_matrix):
    """
    Calculates the time-seris cube for the given raster and design matrix
    """
    # multiply the pixel's set of coefficients by the design mat, then sum
    # all the regression terms
    # NOTE: the coeffs are a numpy array of shape (5,);
    #       the design matrix is a numpy array of shape (365, 5);
    #       the multiplication operates row by row, and elementwise on each
    #       row, returning a new (365, 5) array that, when summed across the
    #       columns (i.e. the regression's linear terms) gives a (365,) vector
    #       of fitted daily values for the pixel
    # NOTE: for now, just adding a length-1 second dimension
    #       in order to match the EOF demo code online
    ts_cube = np.nan * np.ones([desmat.shape[0]]+ [1] + [*rast.shape[1:]],
                               np.float32)
    for i in range(ts_cube.shape[2]):
        for j in range(ts_cube.shape[3]):
            if not np.any(np.isnan(rast[:, i, j])):
                ts_cube[:, 0, i, j] = calc_time_series(rast,
                                                       i, j,
                                                       design_matrix)
    return ts_cube


def calc_time_series(rast, i, j, design_matrix):
    """
    Calculates the time series at pixel i,j, using the coefficients for the
    constant and the sin and cosine terms of the annual and semiannual
    harmonic components. Returns the time series as a numpy array.
    """
    ts = np.sum(rast[:, i, j] * design_matrix, axis=1)
    return ts


rast_file = rio.open('./mosaic.tif')
rast = rast_file.read()

transform = rast_file.get_transform()
dims = rast.shape[1:]

lats, lons = get_lats_lons(transform, dims)

# names of the bands saved into the TFRecord files
band_names = ['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']

# design matrix and time series
desmat = make_design_matrix(thin=True, thin_frac=0.1)

# subset the lons and lats, for memory mgmt
lon_subset_minidx, lon_subset_maxidx = (0, len(lons))
lat_subset_minidx, lat_subset_maxidx = (0, len(lats))

# time series array
ts_cube = calc_time_series_cube(rast[:,
                                     lat_subset_minidx:lat_subset_maxidx,
                                     lon_subset_minidx:lon_subset_maxidx],
                                desmat)

# use weights? (cosine of latitude)
use_weights = False


# DETH: 08-16-21: code adapted from:
    # https://ajdawson.github.io/eofs/latest/examples/nao_standard.html
# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights can be applied before the computation of EOFs.
if use_weights:
    coslat = np.cos(np.deg2rad(lats)).clip(0., 1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]
    solver = Eof(ts_cube, weights=wgts)
else:
    solver = Eof(ts_cube)


clevs = np.linspace(-75, 75, 11)
proj = ccrs.Orthographic(central_longitude=np.mean(lons),
                         central_latitude=np.mean(lats))
# Retrieve the 2 leading EOFs, expressed as the covariance between the leading PC
# time series and the input fitted SIF/NIRvP values
eof1 = solver.eofsAsCovariance(neofs=2)[0,:,:,:]
eof2 = solver.eofsAsCovariance(neofs=2)[1,:,:,:]

for eof_n, title_n in [(eof1, 'EOF1'), (eof2, 'EOF2')]:
    fig = plt.figure()
    ax = plt.axes(projection=proj)
    # Plot the leading EOF expressed as covariance
    ax.set_global()
    ax.coastlines()
    ax.contourf(lons[lon_subset_minidx:lon_subset_maxidx],
                lats[lat_subset_minidx:lat_subset_maxidx],
                eof_n.squeeze(), levels=clevs,
                cmap=plt.cm.coolwarm, transform=ccrs.PlateCarree())
    plt.title('%s expressed as covariance' % title_n, fontsize=16)
    plt.show()


# show variance fractions (i.e. scree plot)
varfrac = solver.varianceFraction()
plt.figure()
eof_num = range(1, len(varfrac)+1)
plt.plot(eof_num, varfrac, linewidth=2)
plt.plot(eof_num, varfrac, linestyle='None', marker="o", color='r', markersize=8)
plt.axhline(0, color='k')
plt.xticks(range(1, len(varfrac)+1))
plt.title('Fraction of the total variance represented by each EOF')
plt.xlabel('EOF #')
plt.ylabel('Variance Fraction')
plt.xlim(1, len(varfrac))
plt.ylim(np.min(varfrac), np.max(varfrac)+0.01)
plt.show()

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
from sklearn.cluster import KMeans
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from copy import deepcopy
from pprint import pprint
import rasterstats as rs
import rioxarray as rxr
import geopandas as gpd
import rasterio as rio
import pandas as pd
import numpy as np
import glob
import json
import time
import re
import os

import matplotlib.pyplot as plt
from shapely.ops import unary_union
from shapely.geometry import Polygon, Point
from scipy.spatial import ConvexHull
from collections import Counter as C



####################
# VARS AND FILEPATHS
####################

# set the neighborhood radii (in km)
neigh_rads = [50, 100, 150]

# get the variables' filepaths
REPO_DIR = ('/home/deth/Desktop/CAL/research/projects/seasonality/'
            'seasonal_asynchrony')
DATA_DIR = ('/home/deth/Desktop/CAL/research/projects/seasonality/'
            'seasonal_asynchrony/data')
FIGS_DIR = ('/home/deth/Desktop/CAL/research/projects/seasonality/'
            'seasonal_asynchrony/res/figs')
TABS_DIR = ('/home/deth/Desktop/CAL/research/projects/seasonality/'
            'seasonal_asynchrony/res/tabs')
EXTERNAL_DATA_DIR = '/media/deth/SLAB/diss/3-phn/final_maps_and_results/'
EXTERNAL_MASK_DATA_DIR = '/media/deth/SLAB/diss/3-phn/GEE_outputs/LSP_masks/'
EXTERNAL_RF_DATA_DIR = '/media/deth/SLAB/diss/3-phn/final_maps_and_results/rf/'
EXTERNAL_INAT_DATA_DIR = '/media/deth/SLAB/diss/3-phn/inat/'
EXTERNAL_FLUX_DATA_DIR = '/media/deth/SLAB/diss/3-phn/flux/'
COEFFS_FILE = os.path.join(EXTERNAL_DATA_DIR, 'NIRv_coeffs.tif')
COEFFS_STRICT_FILE = os.path.join(EXTERNAL_DATA_DIR, 'NIRv_STRICT_coeffs.tif')
EOFS_FILE = os.path.join(EXTERNAL_DATA_DIR,
                         'NIRv_4_EOFs_sqrt_coswts_standts.tif')
EOFS_PREPPED_FILE = os.path.join(EXTERNAL_DATA_DIR,
                         'NIRv_4_EOFs_sqrt_coswts_standts_SCALED_FOLDED_EPSG-8857.tif')
ASYNCH_FILES = {rad: os.path.join(EXTERNAL_DATA_DIR,
                'NIRv_STRICT_asynch_%ikm.tif' % rad) for rad in neigh_rads}
BOUNDS_DIR = os.path.join(DATA_DIR, 'bounds')
ADM0_BOUNDS = os.path.join(BOUNDS_DIR, 'global_adm0.shp')
GLOBAL_ADM1_BOUNDS = os.path.join(BOUNDS_DIR, 'global_adm1.shp')
SELECT_ADM1_BOUNDS = os.path.join(BOUNDS_DIR, 'select_adm1.shp')
BIOCLIM_DIR = os.path.join(DATA_DIR, 'bioclim')
BIOCLIM_INFILEPATHS = glob.glob(os.path.join(BIOCLIM_DIR,"wc2.1_2.5m_bio_*.tif"))
BIOCLIM_INFILEPATHS = sorted(BIOCLIM_INFILEPATHS)

# pattern that occurs just before the file number in each input file's name
PATT_B4_FILENUM = '-OUT-'

# kernel size used by GEE to output the TFRecord files
KERNEL_SIZE = 60

# default missing-data val
DEFAULT_VAL = -9999.0

# names of the bands saved into the TFRecord files
INBANDS = ['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']
OUTBANDS = ['asynch', 'asynch_R2', 'asynch_euc', 'asynch_euc_R2', 'asynch_n']



###############################
# HARMONIC REGRESSION FUNCTIONS
###############################

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



#################################
# ABSTRACT QUANTITATIVE FUNCTIONS
#################################

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



####################
# TEMPORAL FUNCTIONS
####################

def calc_doy_diff(doy1, doy2):
    '''
    calculate the distance, in number of days, between 2 numericaly days of year
    '''
    # get the lesser of the distance back to the earlier day of year or
    # forward to the same day of year next year
    d = sorted([doy1, doy2])
    dist = np.min((d[1]-d[0], d[0] + 365 - d[1]))
    return dist



######################
# GEOSPATIAL FUNCTIONS
######################

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


def make_xarr_bbox_poly(xarr, n_xcoords=4000, n_ycoords=2000):
    xmin, ymin, xmax, ymax = xarr.rio.bounds()
    coords = []
    # draw closely spaced coordinates all along the perimeter
    # (so that projected box, in Equal Earth proj, isn't a simple trapezoid
    # that coarsely excises pixels in the wider equatorial regions of the
    # projection)
    coords.extend([(xmin, n) for n in np.linspace(ymin, ymax, n_ycoords)])
    coords.extend([(n, ymax) for n in np.linspace(xmin, xmax, n_xcoords)])
    coords.extend([(xmax, n) for n in np.linspace(ymax, ymin, n_ycoords)])
    coords.extend([(n, ymin) for n in np.linspace(xmax, xmin, n_xcoords)])
    bbox_poly = Polygon(coords)
    return bbox_poly


def mask_xarr_to_other_xarr_bbox(xarr,
                                 other_xarr,
                                 drop=False,
                                 n_bbox_xcoords=4000,
                                 n_bbox_ycoords=2000,
                                ):
    '''
    mask the first xarray object to the bounding box of the second
    (NOTE: CRS differences will be reconciled automatically!)
    '''
    mask_bbox = make_xarr_bbox_poly(other_xarr,
                                    n_xcoords=n_bbox_xcoords,
                                    n_ycoords=n_bbox_ycoords,
                                   )
    xarr_masked = xarr.rio.clip([mask_bbox],
                                crs=other_xarr.rio.crs,
                                drop=drop,
                               )
    return xarr_masked


def get_raster_info_points(rast_filepath, pts, return_format='vals',
                           standardize=False, fill_nans=False, fill_tol=5,
                           minmax_scale_rast=False):
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

    # min-max scale the raster, if requested
    if minmax_scale_rast:
        rast = (rast - np.min(rast))/(np.max(rast) - np.min(rast))

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


def draw_random_points_within_polygon(n, polygon, verbose=False):
    '''
    faster implementation of drawing random points within polygon,
    borrowed from https://www.matecdev.com/posts/random-points-in-polygon.html
    '''
    minx, miny, maxx, maxy = polygon.bounds
    out_gdfs = []
    total_n_out = 0
    i=0
    while total_n_out < n:
        x = np.random.uniform(minx, maxx, 2*n)
        y = np.random.uniform(miny, maxy, 2*n)
        df = pd.DataFrame()
        df['pts'] = list(zip(x,y))
        df['pts'] = df['pts'].apply(Point)
        gdf_pts = gpd.GeoDataFrame(df, geometry='pts')
        gdf_poly = gpd.GeoDataFrame(index=['target_poly'], geometry=[polygon])
        sjoin = gpd.tools.sjoin(gdf_pts, gdf_poly, predicate='within', how='left')
        pts_in_poly = gdf_pts[sjoin.index_right=='target_poly']
        out_gdfs.append(pts_in_poly.iloc[:n, :])
        total_n_out += len(pts_in_poly.iloc[:n, :])
        i += 1
    if verbose:
        print(f"\tdraw_random_points_within_polygon ran {i} iterations of the loop")
    out_pts = pd.concat(out_gdfs).iloc[:n, :]
    return out_pts


def calc_pw_clim_dist_mat(pts,
                          nodata_val=-3.4e+38,
                          vars_list=None,
                         ):
    """
    Calculates the matrix of pairwise distances between
    any (if vars_list is not None) or all standardized
    BioClim variables for the given nx2 numpy array
    geographic coordinates (with the columns provided in lon,lat order).
    """
    # get only the indicated subset of bioclim vars, if necessary
    if vars_list is None:
        infilepaths = [*BIOCLIM_INFILEPATHS]
    else:
        find_var = lambda s: re.search('bio_\d{1,2}(?=\.tif)', s).group()
        infilepaths = [fp for fp in BIOCLIM_INFILEPATHS if find_var(fp) in vars_list]

    # list to hold arrays of all bioclim vars' vals
    bioclim_vals_arrs = []
    for fn in infilepaths:
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
    assert bioclim.shape[0] == pts.shape[0]
    assert bioclim.shape[1] == len(infilepaths)
    if vars_list is not None:
        assert bioclim.shape[1] == len(vars_list)
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



####################
# PLOTTING FUNCTIONS
####################

def strip_axes_labels_and_ticks(ax):
    '''
    get rid of axis labels, ticks, and title
    '''
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('')


def plot_juris_bounds(ax=None,
                      lev0=True,
                      lev0_color='none',
                      lev0_linecolor='black',
                      lev0_linewidth=0.5,
                      lev0_alpha=0.7,
                      lev0_zorder=2,
                      lev1=True,
                      lev1_which='large',
                      lev1_color='none',
                      lev1_linecolor='black',
                      lev1_linewidth=0.3,
                      lev1_alpha=0.5,
                      lev1_zorder=1,
                      crs=8857,
                      strip_axes=True,
                      reset_axlims=False
                     ):
    '''
    plot jurisdictional boundaries on the given axes,
    including administrative level 0 and/or level 1 boundaries,
    and including either all global level 1 boundaries or
    just select nations with large land areas (Canada, US, Mexico, Brazil,
    Argentina, Kazakhstan, India, Russia, China, Australia)
    '''
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        return_ax = True
    else:
        return_ax = False
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
    if lev0:
        # load country boundaries
        adm0 = gpd.read_file(ADM0_BOUNDS).to_crs(crs)
        adm0.plot(color=lev0_color,
                  edgecolor=lev0_linecolor,
                  linewidth=lev0_linewidth,
                  alpha=lev0_alpha,
                  zorder=lev0_zorder,
                  ax=ax,
                 )
    if lev1:
        assert lev1_which in ['large', 'all']
        lev1_file = {'large': SELECT_ADM1_BOUNDS,
                     'all': GLOBAL_ADM1_BOUNDS}[lev1_which]
        adm1 = gpd.read_file(lev1_file).to_crs(crs)
        adm1.plot(color=lev1_color,
                  edgecolor=lev1_linecolor,
                  linewidth=lev1_linewidth,
                  alpha=lev1_alpha,
                  zorder=lev1_zorder,
                  ax=ax,
                 )
    if strip_axes:
        strip_axes_labels_and_ticks(ax)
    if return_ax:
        return ax
    else:
        if not reset_axlims:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)



def plot_flowerdate_LSP_comparison(flower_obs,
                                   ax_map,
                                   ax_ts,
                                   ax_radar,
                                   colors=np.array(['#2d5098', '#ca1957']), # [blue, red]
                                   plot_crs=8857,
                                   map_xlim=None,
                                   map_ylim=None,
                                   interp_lsp_data=False,
                                   neigh_dist_lsp_fill_tol=2,
                                   radar_alpha=0.5,
                                   radar_width_shrink_factor=0.9,
                                   save_scree_plot=False,
                                   name=None,
                                   tid=None,
                                  ):
    """
    plot a visual comparison between the flowering-doy differences
    across a set of iNat observation locations and the spatial LSP variability
    observed at the locations, using the length of the 'colors' argument
    to determine K (i.e., the number of clusters to fit to the LSP
    time series data using K-means clustering)
    """
    # get the LSP time series at all observations
    pts = flower_obs.get_coordinates().values
    lsp_ts = get_raster_info_points(COEFFS_STRICT_FILE,
                                     pts,
                                     'ts',
                                     standardize=True,
                                     fill_nans=interp_lsp_data,
                                     fill_tol=neigh_dist_lsp_fill_tol,
                                    )
    assert np.all(lsp_ts.shape == (len(flower_obs), 365))
    # TODO: determine and drop those without those LSP data available
    miss = np.sum(pd.isnull(lsp_ts), axis=1) > 0
    pct_lsp_miss = np.mean(miss)
    pct_lsp_miss_msg = (f"{np.round(100*pct_lsp_miss, 2)}% "
                        f"of sites for are missing LSP data")
    print(pct_lsp_miss_msg)
    keep = np.invert(miss)
    flower_obs = flower_obs.loc[keep, :]
    pts = pts[keep, :]
    lsp_ts = lsp_ts[keep, :]
    assert np.all(lsp_ts.shape == (len(flower_obs), 365))
    assert np.all(pts.shape == (len(flower_obs), 2))
    # cluster by genetics and by LSP
    K = len(colors)
    lsp_clusts = KMeans(n_clusters=K).fit(lsp_ts)
    # also save a scree plot for values of K=[1,10], if requested
    if save_scree_plot:
        assert name is not None and tid is not None
        fig_scree, ax_scree = plt.subplots(1,1)
        inertias = [KMeans(n_clusters=i).fit(lsp_ts).inertia_ for i in range(1,
                                                                             11)]
        ax_scree.plot([*range(1, 11)], inertias)
        ax.set_xlabel('K')
        ax.set_ylabel('K-means clustering inertia')
        ax.set_title(f"{name}, TID:{tid}")
        fig_scree.savefig(os.path.join(FIGS_DIR,
                                       'MMRR_res_figs',
            f"TID_{tid}_{name.replace(' ', '_')}_lsp_ts_kmeans_clust_scree.png"),
                          dpi=300,
                         )
    # get the first index of the points at the northernmost
    # location in the dataset, as this will be used to ensure that
    # that point's cluster is blue in each figure (for clearer visual
    # comparison)
    northest_pt_idx = np.where(pts[:, 1] == np.max(pts[:, 1]))[0][0]
    # flip the clustering labels
    # so that the first of the northernmost points is always labeled
    # cluster 0 (because the labels are arbitrary anyhow,
    # but this aligns the colors with Thomé Fig. 5,
    # and aligns them across our subfigs, by putting blue in north)
    if lsp_clusts.labels_[northest_pt_idx] == 1:
        lsp_clusts.labels_ = np.int16(lsp_clusts.labels_ == 0)
    # jitter points by <=0.05 degrees, i.e., ~<=1 of our analysis raster pixels,
    # for visibility of co-ocurring samples
    x_jitters = np.random.uniform(low=-0.05, high=0.05, size=pts.shape[0])
    y_jitters = np.random.uniform(low=-0.05, high=0.05, size=pts.shape[0])
    jitters = np.stack((x_jitters, y_jitters)).T
    assert np.all(pts.shape == jitters.shape)
    pts_for_plot = pts + jitters
    # project the plotting plots as needed
    pts_df = pd.DataFrame({'geometry': [Point(*pts_for_plot[i,:]) for i in
                                        range(pts_for_plot.shape[0])],
                           'idx': [*range(pts_for_plot.shape[0])]
                          })
    pts_gdf = gpd.GeoDataFrame(pts_df, geometry='geometry', crs=4326)
    pts_for_plot = pts_gdf.to_crs(plot_crs).get_coordinates().values
    # map sample points, colored by both the LSP and genetic clusters
    ax_map.scatter(pts_for_plot[:, 0],
                   pts_for_plot[:, 1],
                   s=11,
                   c=colors[lsp_clusts.labels_],
                   edgecolor='black',
                   linewidth=0.5,
                   alpha=0.6,
                   zorder=2,
                  )
    plot_juris_bounds(ax_map,
                      lev1_linewidth=0.05,
                      lev1_alpha=0.6,
                      lev1_zorder=0,
                      lev0_linewidth=0.2,
                      lev0_alpha=1,
                      lev0_zorder=1,
                      crs=plot_crs,
                      strip_axes=True,
                     )
    if map_xlim is not None:
        ax_map.set_xlim(*map_xlim)
    if map_ylim is not None:
        ax_map.set_ylim(*map_ylim)
    # plot LSP fitted time series, colored by both the LSP and genetic clusters
    for i in range(lsp_ts.shape[0]):
        ax_ts.plot(lsp_ts[i, :],
                   color=colors[lsp_clusts.labels_[i]],
                   alpha=0.7,
                   linewidth=0.5,
                  )
    xax_ticks = [0, 90, 180, 271, 364]
    xax_ticklabs = ['Jan', 'Apr', 'Jul', 'Oct', 'Jan']
    ax_ts.set_xticks(xax_ticks)
    ax_ts.set_xticklabels(xax_ticklabs, size=7)
    assert np.round(np.max(lsp_ts), 0) == 2
    assert np.round(np.min(lsp_ts), 0) == -2
    ax_ts.set_yticks([*range(-2, 3)],
                  [str(n) for n in range(-2, 3)],
                  size=7,
                 )
    ax_ts.set_xlim(0, 364)
    ax_ts.set_xlabel('doy', fontdict={'fontsize': 8})
    ax_ts.set_ylabel('scaled LSP', fontdict={'fontsize': 8})
    # plot radar plot of observed flowering dates
    weeks = np.linspace(0, int(7*np.ceil(365/7)), int(np.ceil(365/7)))
    woy = [np.where((weeks[:-1]<=d) * (d<weeks[1:]))[0] for d in flower_obs['doy']]
    assert np.all([len(w) == 1 for w in woy])
    flower_obs.loc[:, 'woy'] = [w[0] for w in woy]
    for lab in np.unique(lsp_clusts.labels_):
        color = colors[lab]
        sub_obs = flower_obs[lsp_clusts.labels_ == lab]
        woy_cts = {i:0 for i in range(52)}
        for i, row in sub_obs.iterrows():
            w = row['woy']
            woy_cts[w] += 1
        for w, ct in woy_cts.items():
            ax_radar.bar(x=w/52*2*np.pi,
                         height=ct,
                         width=2*np.pi/52*radar_width_shrink_factor,
                         alpha=radar_alpha,
                         color=color,
                        )
    # TODO: ADD AXIS LABELS FOR MONTHS


def plot_popgen_LSP_comparison(gen_dist_mat,
                               pts,
                               ax_lspclust_map,
                               ax_lspclust_ts,
                               ax_genclust_map,
                               ax_genclust_ts,
                               colors=np.array(['#2d5098', '#ca1957']), # [blue, red]
                               plot_crs=8857,
                               map_xlim=None,
                               map_ylim=None,
                               interp_lsp_data=False,
                               neigh_dist_lsp_fill_tol=2,
                              ):
    """
    plot a visual comparison between the genetic structure of the given genetic
    distance matrix and the spatial LSP variability observed as the given
    points, using the length of the 'colors' argument to set K
    (i.e., the number of genetic and LSP clusters to fit by K-means clustering)
    """
    # cluster by genetics and by LSP
    K = len(colors)
    gen_clusts = KMeans(n_clusters=K).fit(gen_dist_mat)
    lsp_ts = get_raster_info_points(COEFFS_STRICT_FILE,
                                    pts,
                                    'ts',
                                    standardize=True,
                                    fill_nans=interp_lsp_data,
                                    fill_tol=neigh_dist_lsp_fill_tol,
                                   )
    lsp_clusts = KMeans(n_clusters=K).fit(lsp_ts)
    # get the first index of the points at the northernmost
    # location in the dataset, as this will be used to ensure that
    # that point's cluster is blue in each figure (for clearer visual
    # comparison)
    northest_pt_idx = np.where(pts[:, 1] == np.max(pts[:, 1]))[0][0]
    # flip the clustering labels in both sets of clustering results,
    # so that the first of the northernmost points is always labeled
    # cluster 0 (because the labels are arbitrary anyhow,
    # but this aligns the colors with Thomé Fig. 5,
    # and aligns them across our subfigs, by putting blue in north)
    if gen_clusts.labels_[northest_pt_idx] == 1:
        gen_clusts.labels_ = np.int16(gen_clusts.labels_ == 0)
    if lsp_clusts.labels_[northest_pt_idx] == 1:
        lsp_clusts.labels_ = np.int16(lsp_clusts.labels_ == 0)
    # (jittered by <=0.05 degrees, i.e., ~<=1 of our analysis raster pixels,
    # for visibility of co-ocurring samples)
    x_jitters = np.random.uniform(low=-0.05, high=0.05, size=pts.shape[0])
    y_jitters = np.random.uniform(low=-0.05, high=0.05, size=pts.shape[0])
    jitters = np.stack((x_jitters, y_jitters)).T
    assert np.all(pts.shape == jitters.shape)
    pts_for_plot = pts + jitters
    # project the plotting plots as needed
    pts_df = pd.DataFrame({'geometry': [Point(*pts_for_plot[i,:]) for i in
                                        range(pts_for_plot.shape[0])],
                           'idx': [*range(pts_for_plot.shape[0])]
                          })
    pts_gdf = gpd.GeoDataFrame(pts_df, geometry='geometry', crs=4326)
    pts_for_plot = pts_gdf.to_crs(plot_crs).get_coordinates().values
    # map sample points, colored by both the LSP and genetic clusters
    for ax, clusts in zip([ax_lspclust_map, ax_genclust_map],
                          [lsp_clusts, gen_clusts]):
        ax.scatter(pts_for_plot[:, 0],
                  pts_for_plot[:, 1],
                   s=11,
                   c=colors[clusts.labels_],
                   edgecolor='black',
                   linewidth=0.5,
                   alpha=0.6,
                   zorder=2,
                  )
        plot_juris_bounds(ax,
                         lev1_linewidth=0.05,
                         lev1_alpha=0.6,
                         lev1_zorder=0,
                         lev0_linewidth=0.2,
                         lev0_alpha=1,
                         lev0_zorder=1,
                         crs=plot_crs,
                         strip_axes=True,
                        )
        if map_xlim is not None:
            ax.set_xlim(*map_xlim)
        if map_ylim is not None:
            ax.set_ylim(*map_ylim)
    # plot LSP fitted time series, colored by both the LSP and genetic clusters
    for ax, clusts in zip([ax_lspclust_ts, ax_genclust_ts],
                          [lsp_clusts, gen_clusts],
                         ):
        for i in range(lsp_ts.shape[0]):
            ax.plot(lsp_ts[i, :],
                    color=colors[clusts.labels_[i]],
                    alpha=0.7,
                    linewidth=0.5,
                   )
        xax_ticks = [0, 90, 180, 271, 364]
        xax_ticklabs = ['Jan', 'Apr', 'Jul', 'Oct', 'Jan']
        ax.set_xticks(xax_ticks)
        ax.set_xticklabels(xax_ticklabs, size=7)
        assert np.round(np.max(lsp_ts), 0) == 2
        assert np.round(np.min(lsp_ts), 0) == -2
        ax.set_yticks([*range(-2, 3)],
                      [str(n) for n in range(-2, 3)],
                      size=7,
                     )
        ax.set_xlim(0, 364)
        ax.set_xlabel('day of calendar year', fontdict={'fontsize': 8})
        ax.set_ylabel('normalized LSP value', fontdict={'fontsize': 8})


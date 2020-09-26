#!/bin/python
# calc_asynchrony.py

"""
Reads in the data from all the TFRecord files pertaining to a mixerfile.
Calculates asynchrony for that data and writes out to a matching set of files.

NOTES:
    - TFRecord terminology:
        - 'example' is synonymous with GEE's 'patch'

"""

#--------
# imports
#--------

import matplotlib.pyplot as plt
import geopandas as gpd
from geopy.distance import geodesic
from shapely.geometry import Point
from pprint import pprint
import numpy as np
import tensorflow as tf
import glob
import json
import os
from sklearn.linear_model import LinearRegression


#-----------
# set params
#-----------

# directory where the data and mixerfile live
DATA_DIR = ('/home/drew/Desktop/stuff/berk/research/projects/seasonality/'
              'GEE_output')

# kernel size used by GEE to output the TFRecord files
KERNEL_SIZE = 100

# names of the bands saved into the TFRecord files
BANDS = ['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']

# max distance out to which to find and include neighbors in
# each pixel's asynchrony calculation (in meters)
NEIGH_RAD = 150_000


#-----------------
# define functions
#-----------------

def read_mixer_file(data_dir, return_dims=True):
    """
    Gets the info out of the mixerfile located at data_dir.
    """
    mixerfilepaths = glob.glob(os.path.join(data_dir, '*mixer.json'))
    assert len(mixerfilepaths) == 1, "MORE THAN 1 MIXER FILE FOUND!"
    mixerfilepath = mixerfilepaths[0]

    # read the mixer file
    mixer = json.load(open(mixerfilepath, 'r'))
    return mixer


def get_mixer_info(mixer_content):
    """
    Parses the needed information out of a dict of mixerfile contents,
    altering it as necessary.
    """
    # get patch dimensions, adding the kernel size to correct for error
    # in the GEE TFRecord output
    dims = calc_patch_dimensions(mixer_content, KERNEL_SIZE)
    # get the SRS and its affine projection matrix
    srs = mixer_content['projection']
    crs = srs['crs']
    affine = srs['affine']['doubleMatrix']
    # get the x and y resolutions
    xres, yres = affine[0], affine[4]
    # get the x and y min values (i.e. the center coordinates of
    # the top-left pix)
    # NOTE: need to subtract 50 pixels' worth of res, to account for
    #       fact that mixer file doesn't include the size of the kernel
    #       used to create neighbor-patch overlap
    # NOTE: need to add half pixel's worth of res, to account for
    #       fact that the xmin and ymin are expressed
    #       as upper-left pixel corner, whereas I want to operate
    #       on basis of pixel centers
    xmin = affine[2] - (((KERNEL_SIZE/2) + 0.5)* xres)
    ymin = affine[5] - (((KERNEL_SIZE/2) + 0.5)* yres)
    xmax = xmin + (dims[0] * xres)
    ymax = ymin + (dims[1] * yres)
    # get the number of patches per row and the total number of rows
    # NOTE: FOR NOW, ASSUMING THAT THE MIXER IS SET UP SUCH THAT THE MOSAIC SHOULD
    # BE FILLED ROW BY ROW, LEFT TO RIGHT
    patches_per_row = mixer_content['patchesPerRow']
    tot_patches = mixer_content['totalPatches']
    return dims, crs, xmin, ymin, xres, yres, patches_per_row, tot_patches


def calc_patch_dimensions(mixer_content, kernel_size):
    """
    Calculates and returns the patch dimensions using the input mixerfile info
    and kernel size.
    """
    # Get relevant info from the JSON mixer file
    # NOTE: adding overlap to account for the kernel size,
    # which isn't reflected in the mixer file
    patch_width = mixer_content['patchDimensions'][0] + kernel_size
    patch_height = mixer_content['patchDimensions'][1] + kernel_size
    patch_dimensions = (patch_width, patch_height)
    return patch_dimensions


def parse_tfexample_to_numpy(ex, dims, bands, default_val=-9999.0):
    """
    Parse a TFRecordDataset's Example to a 3D lon x lat x n_bands array,
    then return it.
    """
    arrays = []
    # parse each feature into an identically shaped numpy array
    for beta in bands:
        # get example from TFRecord string
        parsed = tf.train.Example.FromString(ex)
        # pull the feature corresponding to this beta
        feature = parsed.features.feature[beta]
        floats = feature.float_list.value
        arr = np.array(floats).reshape(dims)
        # NOTE: PROBABLY DELETE THE MASKING LINES
        #mask = np.int8(arr == default_val)
        #print(mask)
        #masked = np.ma.array(arr, mask=np.int8(arr == default_val))
        arr[arr == default_val] = np.nan
        arrays.append(arr)
    out = np.stack(arrays)
    return out


def read_tfrecord_file(infilepath, dims, bands):
    """
    Combines all the patches contained in the file located at infilepath
    into a TFRecordDataset. Returns a generator that yields each next example
    (i.e. patch) in that dataset, parsed into an n_bands x lon x lat numpy array.
    """
    # create a TFRecordDataset from the files
    raw_dataset = tf.data.TFRecordDataset(infilepath)
    # turn the dataset into a list of arrays
    for example in raw_dataset.as_numpy_iterator():
        yield parse_tfexample_to_numpy(example, dims, bands)

def get_patch_lons_lats(xmin, ymin, xres, yres, dims, col_j, row_i):
    """
    Takes the overall x and y min values of the TFRecord dataset,
    the x and y resolutions, the patch dimensions, and the column and row
    indices of the current patch. Returns a meshgrid of the cell centers of the
    current patch.
    """
    # calculate the x and y mins of the current patch
    patch_xmin = xmin + (col_j * dims[0] * xres)
    patch_ymin = ymin + (row_i * dims[1] * yres)

    # get lists of xs and ys of all the current patch's pixels
    xs = np.linspace(patch_xmin, xres, dims[0])
    ys = np.linspace(patch_ymin, yres, dims[1])

    # get the meshgrid of those coordinates
    #gridx, gridy = [coords.flatten() for coords in np.meshgrid(xs, ys)]
    gridx, gridy = [coords for coords in np.meshgrid(xs, ys)]
    return gridx, gridy


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


def calc_time_series(patch, i, j, regression_vals):
    """
    Calculates the time series at pixel i,j, using the coefficients for the
    constant and the sin and cosine terms of the annual and semiannual
    harmonic components. Returns the time series as a numpy array.
    """
    # multiply the pixel's set of coefficients by the regression vals, then sum
    # all the regression terms
    # NOTE: the coeffs are a numpy array of shape (5,);
    #       the regression values are a numpy array of shape (365, 5)
    #       the multiplication operates row by row, and elementwise on each
    #       row, returning a new (365, 5) array that, when summed across the
    #       columns (i.e. the regression's linear terms) gives a (365,) vector
    #       of fitted daily values for the pixel
    ts = np.sum(patch[:, i, j] * regres_vals, axis=1)
    return ts


def linear_regression(y, x):
    """
    Calculates the simple linear regression of y on x.

    Returns the regression as a model object.
    """
    # Reshape the y variable
    reshape_y = np.array(y).reshape(-1,1)
    # Make a model
    model = LinearRegression().fit(reshape_y, x)
    return model 


def calc_distance(lon1, lat1, lon2, lat2):
    """
    Calculates the real-world distance (in meters)
    between the two points at lon1,lat1 and lon2,lat2.

    Returns the distance as a float.
    """

    pass


def find_neighborhood_pixels(i, j, ys, xs, rast, yres, xres, crs,
                             neigh_rad=NEIGH_RAD, min_dist_deg_lat=111_000):
    """
    Finds all pixel-center coordinates within neigh_dist (in meters) of
    the input lat,long coordinates.

    To do this, first it uses the minimum distance in meters of a degree
    longitude within our dataset (based on the maximum absolute value of a
    latitudinal point included in the dataset; here set to 67 degrees
    latitude (a bit of an overshoot),
    which based upon https://en.wikipedia.org/wiki/Decimal_degrees#Precision
    gives us 43,496 m distance for a degree longitude) to subset a box of
    pixels of maximum necessary size surrounding the focal pixel.
    Then it casts the pixels as shapely Points, then uses
    Geopandas' distance method to grab all pixels within the subsetted box
    whose centers are <= neigh_rad of the focal cell's center. 

    Returns the neighbors' coordinates and distances from the focal pixel's
    center as a 3-tuple (lats, lons, dists) of lists of floats.
    """
    # get the minimum distance of one degree longitude at the maximum
    # absolute-value latitude in the input ys
    max_abs_lat = np.max(np.abs(ys))
    min_dist_deg_lon = geodesic((0, max_abs_lat), (1, max_abs_lat)).m
    # get number of x cells equiv to min_deg_lon_dist, to use as number of
    # neighbors surrounding the focal pixel in each direction that should be
    # included in the subsetted neighborhood raster
    # (cell_tot = m_tot / (m / deg) / (deg / cell)
    print(i, j, ys, xs, rast, yres, xres, crs)
    neigh_rad_y = int(np.ceil(np.abs(neigh_rad/min_dist_deg_lat/yres)))
    neigh_rad_x = int(np.ceil(np.abs(neigh_rad/min_dist_deg_lon/xres)))
    print(neigh_rad_y, neigh_rad_x)

    # use that number to get all the lats and lons from among which neighbors
    # will be chosen
    valid_ys = ys[i-neigh_rad_y:i+neigh_rad_y+1,
                  j-neigh_rad_x:j+neigh_rad_x+1].flatten()
    valid_xs = xs[i-neigh_rad_y:i+neigh_rad_y+1,
                  j-neigh_rad_x:j+neigh_rad_x+1].flatten()
    print(valid_xs, valid_ys)
    
    # save the focal pixel's centerpoints, and a list of the other
    # centerpoints
    foc = (ys[i,j], xs[i,j])
    pts = [(valid_ys[i], valid_xs[i]) for i in range(len(valid_xs))]
    print(foc)
    print(pts)

    # calculate the geodesic distances
    # NOTE: per the Geopy docs
    # (https://geopy.readthedocs.io/en/latest/#module-geopy.distance)
    # I'm assuming the WGS84 ellipsoid model, and ignoring altitude differences
    # between points; seems perfectly fine as far as I'm concerned
    dists = [geodesic(foc, pt).m for pt in pts]
    dists = [dist for dist in dists if dist <= neigh_rad]
    return dists

    """
    # make lat,lon pairs into point geometries
    pts = [Point(coords) for coords in zip(valid_xs.flatten(),
                                           valid_ys.flatten())]
    gdf = gpd.GeoDataFrame(geometry=pts, crs=crs)
    print(gdf)

    # make a second data frame of the focal point repeated,
    # to use with Geopandas' distance method to get distances of
    # all points from focal point
    buff_gdf = gpd.GeoDataFrame(geometry=[Point(xs[i, j], ys[i, j]) for _ in
                                          range(gdf.shape[0])], crs=crs)

    # calculate the distances
    gdf['dists'] = buff_gdf.distance(gdf)
    return gdf
    """
     

def get_infile_outfile_paths(data_dir):
    """
    Uses the data-directory path provided to get and return lists
    of the input and output files' paths.
    """
    # set up IO paths
    infilepaths = glob.glob(os.path.join(data_dir, '*.tfrecord'))
    outfilepaths = [fp.split('.')[0]+'OUT.tfrecord' for fp in infilepaths]
    return infilepaths, outfilepaths


def write_data(rast, filepath):
    """
    Writes the output raster to disk using the given filepath.
    """
    pass



#------------------------------------------------------------------------------
# calculate asynchrony
#---------------------

#----------------
# set up IO paths
#----------------

infilepaths, outfilepaths = get_infile_outfile_paths(DATA_DIR)

#----------------------
# get the design matrix
#----------------------

design_mat = make_design_matrix()

#---------------
# get mixer info
#---------------

mix = read_mixer_file(DATA_DIR)
(dims, crs, xmin, ymin, xres, yres, patches_per_row,
                                    tot_patches) = get_mixer_info(mix)

#--------------
# load the data
#--------------

# loop over the Tfrecord files, keeping track of the patch number and its
# corresponding row and column indices in the overall patch mosaic
patch_n = 0
row_i = 0
col_j = 0
for infilepath in infilepaths:
    # read the file's data in as a set of examples (i.e. patches)
    patches = read_tfrecord_file(infilepath, dims, BANDS)

    # loop over the patches within the current file
    # NOTE: each patch is of shape (n_bands, lat, lon)
    for patch in patches:

        # get the lons and lats of the current example's patch
        xs, ys = get_patch_lons_lats(xmin, ymin, xres, yres, dims, col_j, row_i)

        # create the output asynch array, to be filled in (starts as all NaNs)
        asynch = np.nan * patch[0,:,:]
        
        
        assert True == False
        #----------------
        # run calculation
        #----------------

        # loop over pixels
        for i in patch.shape[1]:
            for j in patch.shape[2]:
        
                # create lists of R2 and dist values
                R2s = []
                dists = []
        
                # calculate the focal pixel's time series
                ts_foc = calc_time_series(patch, i, j, regres_vals)
       
                # TODO: START HERE!
                # NOTES:
                    # - seem to be pulling waaaayy too many candidate pixels
                    # using the approach I'm using now --> debug!
                    # - consider just creating a gpd df, then reprojecting to
                    # an equal-distance CRS and calculating my dists there

                # find the neighborhood of pixels for the focal pixel
                neigh_lats, neigh_lons = find_neighborhood_pixels(lat, lon)
        
                # loop over neighbors
                for ni, nlat in enumerate(neigh_lats):
                    for nj, nlon in enumerate(neigh_lons):
        
                        # get the neighbor's time series
                        ts_neigh = calc_time_series(patch, ni, nj, regres_vals)
        
                        # correlate the two time series and extract the R2
                        R2 = linear_regression(ts_neigh, ts_foc)['R2']
                        R2s.append(R2)
        
                        # get the distance between the focal and neighbor pixels
                        dist = calc_distance(lat, lon, nlat, nlon)
                        dists.append(dist)
                
                # get the slope of the overall regression
                asynch_val = np.abs(linear_regression(R2s, dists)['slope'])
                output[i, j] = asynch_val

        #increment counters
        patch_n += 1
        if col_j == patch_per_row - 1:
            row_i += 1
            col_j = 0
        else:
            col_j += 1

#------------
# data output
#------------

write_data(asynch, outfilepath)


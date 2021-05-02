#!/usr/bin/env python
# coding: utf-8

# In[290]:


import geopandas as gpd
import tensorflow as tf
import numpy as np
import glob
import json
import time
import os
from osgeo import gdal
import random
import numpy as np
from scipy.spatial import KDTree
import math
import itertools
import rasterio as rio


# In[61]:


# BioClim variable
BIO_DATA_DIR = ('C:\\Users\\thaon\\Documents\\Asynchrony\\bioclim')
bio_infilepaths = glob.glob(os.path.join(BIO_DATA_DIR,"*.tif"))
bio_infilepaths = sorted(bio_infilepaths)


# In[284]:


rangeX = (-100, -50) # longtitude
rangeY = (-50, 20) # Latitude
qty = 15  # Enter in the number greater than random points you need

#Generate random x,y coordinates
randPoints = [] # a list to store coordinates
while len(randPoints) < qty:
    x = random.randrange(*rangeX) 
    y = random.randrange(*rangeY) 
    randPoints.append((x,y))
#Extract value at the coordinate
a = [] #a list to store bioclimatic variables
for point in randPoints:
    a += [point]
    for path in bio_infilepaths:
        raster = gdal.Open(path)
        cols = raster.RasterXSize
        rows = raster.RasterYSize
        transform = raster.GetGeoTransform()
        xOrigin = transform[0]
        yOrigin = transform[3]
        pixelWidth = transform[1]
        pixelHeight = -transform[5]   
        r = raster.GetRasterBand(1).ReadAsArray(0,0,cols,rows)
        col = int((point[0] - xOrigin) / pixelWidth)
        row = int((yOrigin - point[1] ) / pixelHeight)
        a +=[r[row][col]]
randPoints #random generated coordinates


# In[285]:


# Functions return list of climate Euclidean distance between pairwise pixel
def filter_out_ocean_pixel(biolist,coords):
    z = np.float32(-3.4e+38) #value of ocean region
    bl = np.array_split(biolist, qty)
    d = np.where(np.all(np.delete(bl,0,1) == np.repeat(z,19),axis=1)) #index of array with pixel over the ocean
    index = np.array(d).tolist() #index to delete
    flat_list = [item for sublist in index for item in sublist]
    for i in sorted(flat_list, reverse=True): #delete array with ocean pixel
        del bl[i]
    ncoords = coords
    for j in sorted(flat_list, reverse=True):
        del ncoords[j]
    return bl, ncoords
def calculate_euc_clim(biolist):
    pw = list(itertools.combinations(list(range(0,len(biolist))),2)) #generate list of pairwise combinations
    ed = []
    for i in pw:
        dist = np.sqrt(np.sum((np.delete(biolist[i[0]],0) - np.delete(biolist[i[1]],0))**2))
        ed += [dist]
    return ed


# In[286]:


# Results
b, coords = filter_out_ocean_pixel(a, randPoints)
ed = calculate_euc_clim(b)


# In[296]:


def get_seasonality_info_points(coeffs_rast_filepath, pts, dists=True):
    """
    takes a raster of coeffs from the harmonic seasonality regression
    and an nx2 np.array of n points' x and y coordinates,
    returns either the fitted time series at those points
    (if 'dists' == False), or a matrix of pairwise seasonal distances
    at those points (if 'dists' == True; default)
    """
    # read in the raster file
    f = rio.open(coeffs_rast_filepath)
    rast_coeffs = f.read()
    print(rast_coeffs.shape)

    # get the cells' max lon- and lat-bound values
    cell_max_lons, cell_max_lats = get_cell_lonlat_bounds(f)
    print('lons', cell_max_lons)
    print('lats', cell_max_lats)

    # get the regression's design matrix
    design_mat = make_design_matrix()

    # matrix-multiply design_mat x coeffs to get fitted time series for all pts
    ts_mat = np.zeros((pts.shape[0], 365)) * np.nan

    for row_i in range(pts.shape[0]):
        # get the array coords of the cell the point falls in
        pt_cell_i, pt_cell_j = get_cell_coords_for_pt(lon=pts[row_i,0],
                                                      lat=pts[row_i,1],
                                                    cell_max_lons=cell_max_lons,
                                                    cell_max_lats=cell_max_lats)
        print('pt_coeffs', rast_coeffs[:, pt_cell_i, pt_cell_j])
        ts = np.sum(rast_coeffs[:, pt_cell_i, pt_cell_j] * design_mat, axis=1)
        ts_mat[row_i, :] = ts

    if dists:
        # calculate pairwise Euc dist matrix
        pw_dist_mat = np.zeros((pts.shape[0], pts.shape[0])) * np.nan
        for row_i in range(pw_dist_mat.shape[0]):
            for col_j in range(pw_dist_mat.shape[1]):
                if row_i == col_j:
                    dist = 0
                else:
                    dist = calc_euc_dist(ts_mat[row_i,:], ts_mat[col_j, :])
                pw_dist_mat[row_i, col_j] = dist
        return pw_dist_mat
    else:
        # return the time series matrix, if pairwise distances not requested
        return ts_mat
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
coeffs_rast_filepath = 'C:\\Users\\thaon\\Documents\\Asynchrony\\SIF_coeffs.tif'
get_seasonality_info_points(coeffs_rast_filepath, coords, dists=True) #how should i store my points?


# In[309]:


n = np.array(coords)
ns = np.split(n,9)
n.shape(0)


# In[129]:


if os.path.abspath('.').split('\\')[1] == 'Users':
    DATA_DIR = ('C:\\Users\\thaon\\Documents\\Asynchrony\\TFP')
else:
    DATA_DIR = ('C:\\Users\\thaon\\Documents\\Asynchrony\\TFP\\')
# pattern that occurs just before the file number in each input file's name
#PATT_B4_FILENUM = 'SIF-'
PATT_B4_FILENUM = 'Amer-'
#PATT_B4_FILENUM = 'NIRvP-'

# kernel size used by GEE to output the TFRecord files
KERNEL_SIZE = 358
#KERNEL_SIZE = 60

# whether or not to trim the half-kernel-width margin before output
TRIM_MARGIN = True

# default missing-data val
DEFAULT_VAL = -9999.0

# names of the bands saved into the TFRecord files
INBANDS = ['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']
OUTBANDS = ['asynch', 'asynch_R2', 'asynch_euc', 'asynch_euc_R2', 'asynch_n']

# max distance out to which to find and include neighbors in
# each pixel's asynchrony calculation (in meters)
NEIGH_RAD = 150_000

# stdout options
VERBOSE = True
TIMEIT = True


# In[130]:


#Function copied from asychrony calc
def get_infile_outfile_paths(data_dir):
    """
    Uses the data-directory path provided to get and return lists
    of the input and output files' paths.
    """
    # set up IO paths
    infilepaths = glob.glob(os.path.join(data_dir, '*.tfrecord'))
    # order the infilepaths
    infilepaths = sorted(infilepaths)
    #make sure that no previously output files are included
    infilepaths = [f for f in infilepaths if not '_OUT' in f]
    # get the corresponding outfilepaths
    outfilepaths = [fp.split('.')[0]+'_OUT.tfrecord' for fp in infilepaths]
    return infilepaths, outfilepaths


def read_mixer_file(data_dir):
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
    patches_per_row = mixer_content['patchesPerRow']
    tot_patches = mixer_content['totalPatches']
    xmin = affine[2] - (((KERNEL_SIZE/2) + 0.5)* xres)
    ymin = affine[5] - (((KERNEL_SIZE/2) + 0.5)* yres)
    xmax = xmin + (dims[0] * xres) * patches_per_row + (KERNEL_SIZE * xres)
    ymax = ymin + (dims[1] * yres) * patches_per_row + (KERNEL_SIZE * yres)
    # get the number of patches per row and the total number of rows
    # NOTE: FOR NOW, ASSUMING THAT THE MIXER IS SET UP SUCH THAT THE MOSAIC SHOULD
    # BE FILLED ROW BY ROW, LEFT TO RIGHT
    
    return dims, crs, xmin, ymin, xmax, ymax, xres, yres, patches_per_row, tot_patches


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


def parse_tfexample_to_numpy(ex, dims, bands, default_val=DEFAULT_VAL):
    """
    Parse a TFRecordDataset's Example to a 3D lon x lat x n_bands array,
    then return it.
    """
    arrays = []
    # get example from TFRecord string
    parsed = tf.train.Example.FromString(ex)
    # coerce each feature into an identically shaped numpy array
    for beta in bands:
        # pull the feature corresponding to this beta
        feature = parsed.features.feature[beta]
        floats = feature.float_list.value
        arr = np.array(floats).reshape(dims)
        arr[arr == default_val] = np.nan
        arrays.append(arr)
    out = np.stack(arrays)
    return out


def read_tfrecord_file(infilepath, dims, bands):
    """
    Combines all the patches contained in the file located at infilepath
    into a TFRecordDataset. Returns a generator that yields each next example
    i.e. patch) in that dataset, parsed into an
    n_bands x lon x lat numpy array.
    """
    # create a TFRecordDataset from the files
    raw_dataset = tf.data.TFRecordDataset(infilepath)
    # turn the dataset into a list of arrays
    for example in raw_dataset.as_numpy_iterator():
        yield parse_tfexample_to_numpy(example, dims, bands)
def serialize_tfrecord_example(patch, bands):
    """
    Serialize a patch (i.e. raster), for writing into a TFRecord file.
    Taken from:
        https://towardsdatascience.com/
            working-with-tfrecords-and-tf-train-example-36d111b3ff4d
    """
    # make the feature dict, with one item for each band label
    # (just like the 'sin_1', 'cos_1', ... bands in the input files)
    feature = {band: tf.train.Feature(float_list=tf.train.FloatList(
               value=patch[i, :, :].flatten())) for i, band in enumerate(bands)}
    #  create a Features message (i.e. protocol buffer) using tf.train.Example
    example_prot = tf.train.Example(features=tf.train.Features(feature=feature))
    return example_prot.SerializeToString()
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


# In[131]:


## Not finished with this code, waiting for timeseries function
mix = read_mixer_file(DATA_DIR)
(DIMS, CRS, XMIN, YMIN, XMAX, YMAX, XRES, YRES, patches_per_row, tot_patches) = get_mixer_info(mix)
infilepaths, outfilepaths = get_infile_outfile_paths(DATA_DIR)
inpatches = read_tfrecord_file(infilepaths,DIMS,INBANDS)
design_mat=make_design_matrix()
#for inpatch in inpatches:
#    calc_time_series(inpatch,-106.5,21.5,design_mat)
XMIN, YMIN, XMAX, YMAX


# In[ ]:





# In[ ]:





# In[ ]:





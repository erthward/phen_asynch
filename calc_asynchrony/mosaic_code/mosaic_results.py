#!/bin/python
# calc_asynchrony.py

"""
Reads in the results calculated by Julia on Savio, mosaics them together,
then saves the global raster to a new file.

NOTE: This is only necessary because the GEE bug won't read my results files
      correctly...

"""


#--------
# imports
#--------

#from sklearn.linear_model import LinearRegression
from shapely.geometry import Point
#from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from pprint import pprint
#import geopandas as gpd
import tensorflow as tf
import rasterio as rio
#from rasterio.profiles import DefaultGTiffProfile
import numpy as np
import glob
import json
import time
import sys
import re
import os


#-----------
# set params
#-----------

# get the files directory
DATA_DIR = sys.argv[1]
#DATA_DIR = ('/home/deth/Desktop/UCB/research/projects'
#            '/seasonality/GEE_output/NA_agg/')

# output filepath
OUTPUT_FILEPATH = sys.argv[2]
#OUTPUT_FILEPATH = 'test_mosaic.tif'

# whether this represents unprocessed GEE output (coeff bands)
# or processed output (asynch metrics)
DATA_TYPE = sys.argv[3].lower()
assert DATA_TYPE in ['c', 'a'], ('the third argument must either be "c" '
                                 '(for coefficients) or "a" (for '
                                 '"asynchrony").')

# pattern that occurs just before the file number in each file's number
PATT_B4_FILENUM = '-OUT-'

# kernel size used by GEE to output the TFRecord files
KERNEL_SIZE = 60
HKW = int(KERNEL_SIZE/2)

# whether or not to trim the half-kernel-width margin before mosaicking
# NOTE: don't need to, because the margins were trimmed before
#       the files were written out
TRIM_MARGIN = True

# default missing-data val
NODATA_VAL = -9999.0

# names of the bands saved into the TFRecord files
if DATA_TYPE == 'a':
    BANDS = ['asynch', 'asynch_R2', 'asynch_euc', 'asynch_euc_R2', 'asynch_n']
else:
    BANDS = ['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']


#-----------------
# define functions
#-----------------

def get_filepaths(data_dir):
    """
    Uses the data-directory path provided to get
    and return a lists of the file paths
    """
    # set up IO paths
    filepaths = glob.glob(os.path.join(data_dir, '*.tfrecord'))
    # order the infilepaths
    filepaths = sorted(filepaths)
    # make sure all files are output files
    filepaths = [f for f in filepaths if 'OUT' in f]
    return filepaths


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
    xmin = affine[2] - (((KERNEL_SIZE/2) + 0.5)* xres)
    ymin = affine[5] - (((KERNEL_SIZE/2) + 0.5)* yres)
    xmax = xmin + (dims[0] * xres)
    ymax = ymin + (dims[1] * yres)
    # get the number of patches per row and the total number of rows
    # NOTE: FOR NOW, ASSUMING THAT THE MIXER IS SET UP
    # SUCH THAT THE MOSAIC SHOULD BE FILLED ROW BY ROW, LEFT TO RIGHT
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


def parse_tfexample_to_numpy(ex, dims, bands, default_val=NODATA_VAL):
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


def get_row_col_patch_ns_allfiles(data_dir, patt_b4_filenum):
    """
    Return an output dict containing the row, column, and patch numbers,
    and outfile paths (as dict values, organized as subdicts)
    for all files (keys).
    """
    print('building FILES_DICT\n')
    # set the starting row, column, and patch counters
    row_i = 0
    col_j = 0
    patch_n = 0

    # get the mixer file info
    mix = read_mixer_file(DATA_DIR)
    (dims, crs, xmin, ymin, xres, yres,
     patches_per_row, tot_patches) = get_mixer_info(mix)

    # get all the input and output file paths
    filepaths = get_filepaths(DATA_DIR)

    # assert that both lists are sorted in ascending numerical order
    # NOTE: if this is not true then my method for tracking the row, col, and
    # patch numbers will break!
    filenums = np.array([int(f.split(patt_b4_filenum)[1].split(
                                         '.tfrec')[0]) for f in filepaths])
    filenums_plus1 = np.array(range(1, len(filepaths)+1))
    assert np.all((filenums_plus1 - filenums) == 1), ("Filepaths do not "
                                                         "appear to be in "
                                                         "numerical order. "
                                                         "\n\t%s") % str(
                                                                    filepaths)

    # make the output dict
    files_dict = {}

    # loop over the input files, get the requisite info (row, col, and patch
    # ns; output file paths), and store in the dict
    for file_i, filepath in enumerate(filepaths):
        print('getting FILES_DICT info for file %s\n' % filepath)
        # create the subdict 
        file_dict = {}
        # stash the outfile path
        file_dict['row_is'] = []
        file_dict['col_js'] = []
        file_dict['patch_ns'] = []

        # read the file into a TFRecordDataset
        dataset = tf.data.TFRecordDataset(filepath)

        # loop over the patches in the infile, store the row, col, and patch
        # numbers, then correctly increment them
        # NOTE: an example (TFRecord jargon) is the same as a patch (GEE jargon)
        for example in dataset:

            # store nums
            file_dict['row_is'].append(row_i)
            file_dict['col_js'].append(col_j)
            file_dict['patch_ns'].append(patch_n)

            #increment counters
            patch_n += 1
            if col_j == patches_per_row - 1:
                row_i += 1
                col_j = 0
            else:
                col_j += 1

        # add this file's file_dict to the files_dict
        files_dict[filepath] = file_dict

    return files_dict


def make_empty_output_array(patches_per_row, tot_patches, dims, bands):
    ncols = patches_per_row * dims[1]
    nrows = int(tot_patches / patches_per_row) * dims[0]
    output = np.ones((len(bands), nrows, ncols)) * np.nan
    return output


def get_patch_insert_indices(row_i, col_j, dims):
    i_start = row_i * dims[0]
    j_start = col_j * dims[1]
    i_stop = i_start + dims[0]
    j_stop = j_start + dims[1]
    return((i_start, i_stop), (j_start, j_stop))


def write_geotiff(output_filepath, output_arr):
    # set up raster metadata
    meta = {'driver': 'GTiff',
            # NOTE: I think it makes sense to stick with float32,
            # since that's the dtype output in the TFRecord files by GEE...
            'dtype': 'float32',
            'nodata': NODATA_VAL,
            'width': output_arr.shape[2],
            'height': output_arr.shape[1],
            'count': output_arr.shape[0],
            'crs': rio.crs.CRS.from_epsg(int(CRS.split(':')[1])),
            'transform': rio.transform.Affine(*mix['projection']['affine'][
                                                            'doubleMatrix']),
            'tiled': False, ## ?
            'interleave': 'band' ## ?
           }
    # make output file profile
    #prof = DefaultGTiffProfile(count=len(BANDS))
    # correct the profile's height and width, nodata val, & number of bands
    #prof.update(blockxsize = output_arr.shape[1],
    #            blockysize = output_arr.shape[0],
    #            nodata = NODATA_VAL,
    #            count=len(BANDS))

    # open the file connection, unpacking the profile
    with rio.open(output_filepath, 'w', **meta) as f:
        # Write an array as a raster band to a new 8-bit file
        f.write(output_arr.astype(rio.float32))


#------------------------------------------------
# get mixer info and set up output data structure
#------------------------------------------------

mix = read_mixer_file(DATA_DIR)
(DIMS, CRS, XMIN, YMIN, XRES, YRES, patches_per_row,
                                    tot_patches) = get_mixer_info(mix)

OUTPUT = make_empty_output_array(patches_per_row, tot_patches, DIMS, BANDS)

#------------------------------------
# get files' row, col, and patch info
#------------------------------------

FILES_DICT = get_row_col_patch_ns_allfiles(DATA_DIR, PATT_B4_FILENUM)

#--------------------------------------------------------
# loop over all files, patches and write them into output
#--------------------------------------------------------

for filepath, patchinfo in FILES_DICT.items():
    print("INSERTING DATA FOR %s\n" % filepath)
    #print(patchinfo)
    # grab patch info
    row_is = patchinfo['row_is']
    col_js = patchinfo['col_js']
    patch_ns = patchinfo['patch_ns']

    # read the file's data into an iterator
    patches = read_tfrecord_file(filepath, DIMS, BANDS)

    # loop over and grab data
    for patch_n, patch in enumerate(patches):
        # get the OUTPUT indices for this patch's data
        i_inds, j_inds = get_patch_insert_indices(row_is[patch_n],
                                                  col_js[patch_n], DIMS)

        # if these are unprocessed GEE output (i.e. coefficients),
        # then trim the margin
        patch = patch[HKW:-HKW, HKW:-HKW]

        # insert the data
        OUTPUT[:, i_inds[0]:i_inds[1], j_inds[0]:j_inds[1]] = patch


#--------------------
# write out to raster
#--------------------

print('\nWRITING RESULTS TO FILE...\n')
write_geotiff(OUTPUT_FILEPATH, OUTPUT)

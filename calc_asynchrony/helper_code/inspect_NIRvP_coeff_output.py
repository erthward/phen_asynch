#!/bin/python
# asynch_fns.py

"""
Reads in the data from all the TFRecord files pertaining to a mixerfile.
Calculates asynchrony for that data and writes out to a matching set of files.

NOTES:
    - TFRecord terminology:
        Each TFRecord file is composed of multiple 'examples'
        (GEE equiv: 'patches'), each of which is composed of multiple
        'features' (GEE: 'bands'), each of which has both
        a 'key' (GEE: 'band name', e.g. 'sin_1') and a value
        (GEE's actual raster of image data; in TFRecord they are stored
        and reading in as 'float_lists' composed of numerous 'value: #######'
        pairs, where ##### is the actual numerical value in a cell).

        In other words, the TFRecord format looks like:

            file_n.tfrecords:

                example 1:

                    feature {
                        key: "constant"
                        value {
                            float_list {
                                value: -9999.0
                                value: -9999.0
                                value: -9999.0
                                .
                                .
                                .
                                }
                            }
                        key: "sin_1"
                        ...
                        key: "sin_2"
                        ...
                        ...
                        }

                example 2:

                    feature {
                        ...
                        }
                ...

                example last:

                    feature {
                        ...
                        }


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
import geopandas as gpd
import tensorflow as tf
import numpy as np
import glob
import json
import time
import os


#-----------
# set params
#-----------

region = 'VA' # 'CA', 'SA', or 'VA'

# directory where the data and mixerfile live
if os.path.abspath('.').split('/')[1] == 'home':
    DATA_DIR = ('/home/deth/Desktop/stuff/berk/research/projects/seasonality/'
                'GEE_output')
else:
    DATA_DIR = '/global/home/users/drewhart/seasonality/GEE_output/SIF/'

# pattern that occurs just before the file number in each input file's name
#PATT_B4_FILENUM = 'SIF-'
#PATT_B4_FILENUM = 'Amer-'
PATT_B4_FILENUM = 'NIRvP-'

# kernel size used by GEE to output the TFRecord files
KERNEL_SIZE = 714

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


#-----------------
# define functions
#-----------------


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


def write_tfrecord_file(patches, filepath, bands):
    """
    Writes the output patches (i.e. rasters) to disk as a TFRecord file,
    using the given filepath.

    Taken from:
        https://towardsdatascience.com/
            working-with-tfrecords-and-tf-train-example-36d111b3ff4d
    """
    #set all NaNs back to the missing-data default val
    for patch in patches:
        patch[np.isnan(patch)] = DEFAULT_VAL
    with tf.io.TFRecordWriter(filepath) as writer:
        for patch in patches:
            # serialize to a TF example
            example = serialize_tfrecord_example(patch, bands)
            # write to disk
            writer.write(example)


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

(dims, crs, xmin, ymin, xres, yres, patches_per_row,
 tot_patches) = get_mixer_info(read_mixer_file('../../GEE_output/%s/' % region))

tfrec = read_tfrecord_file(('../../GEE_output/%s/%s_coeffs_NIRvP-00000.'
                            'tfrecord') % (region, region), dims, INBANDS)

patches = [p for p in tfrec]

plt.imshow(patches[0][1,:,:])
plt.show()

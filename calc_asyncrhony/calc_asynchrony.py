#--------
# imports
#--------

import matplotlib.pyplot as plt
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

MAX_DIST = 300000  # meters

DATA_DIR = ('/home/drew/Desktop/stuff/berk/research/projects/seasonality/'
              'GEE_output')

BANDS = ['beta_sin1', 'beta_cos1', 'beta_sin2', 'beta_cos2']
#BANDS = ['constant', 'beta_sin1', 'beta_cos1', 'beta_sin2', 'beta_cos2']


#-----------------
# define functions
#-----------------

def calc_time_series(i, j, bands, harmonic_terms):
    """
    Calculates the time series at pixel i,j, using the coefficients for the
    sin and cosine terms of the annual and semiannual harmonic components.

    Returns the time series as a list.
    """
    beta_sin1, beta_cos1, beta_sin2, beta_cos2 = bands
    sin1, cos1, sin2, cos2 = harmonic_terms
    ts = (
    #ts = (constant[i, j] +
          beta_sin1[i, j] * sin1 +
          beta_cos1[i, j] * cos1 +
          beta_sin2[i, j] * sin2 +
          beta_cos2[i, j] * cos2)
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


def find_neighborhood_pixels(lat, lon, max_dist=MAX_DIST):
    """
    Finds all pixel-center coordinates within max_dist (in meters) of
    the input lat,long coordinates.

    Returns the neighbors' coordinates as a 2-tuple (lats, lons) of lists of
    floats.
    """
    pass


def write_data(rast, filepath):
    """
    Writes the output raster to disk using the given filepath.
    """
    pass


#------------------
# handle input data
#------------------

# handle filepaths
mixerfilepaths = glob.glob(os.path.join(DATA_DIR, '*mixer.json'))
assert len(mixerfilepaths) == 1, "MORE THAN 1 MIXER FILE FOUND!"
mixerfilepath = mixerfilepaths[0]

infilepaths = glob.glob(os.path.join(DATA_DIR, '*.tfrecord'))

outfilepaths = [fp.split('.')[0]+'OUT.tfrecord' for fp in infilepaths]


# the following data-parsing stuff is with help from:
#        https://colab.research.google.com/github/google/
#        earthengine-api/blob/master/python/examples/ipynb/
#        TF_demo1_keras.ipynb#scrollTo=x2Q0g3fBj2kD)
# and:
#        https://www.tensorflow.org/tutorials/load_data/tfrecord

# read the mixer file
mixer = json.load(open(mixerfilepath, 'r'))

# Get relevant info from the JSON mixer file
# NOTE: add 100 to account for the kernel size,
# which isn't reflected in the mixer file
patch_width = mixer['patchDimensions'][0] + 100
patch_height = mixer['patchDimensions'][1] + 100
patches = mixer['totalPatches']
patch_dimensions = [patch_width, patch_height]
patch_dimensions_flat = [patch_width * patch_height, 1]
#patch_dimensions_flat = [1]

# build the feature-description, for parsing
# NOTE: the tensors are in the shape of a patch, one patch for each band
image_columns = [
      tf.io.FixedLenFeature(shape=patch_dimensions,
#                            dtype=tf.float32,
#                            default_value=-9999.0)
                            dtype=tf.float32)
        for k in BANDS
]
bands_dict = dict(zip(BANDS, image_columns))

# Note that you can make one dataset from many files by specifying a list
dataset = tf.data.TFRecordDataset(infilepath, compression_type='')

# data-parsing function from:
def parse_image(example_proto):
      return tf.io.parse_single_example(tf.io.parse_tensor(example_proto,
                                                           'float32')[0],
                                        bands_dict)

# Parse the data into tensors, one long tensor per patch
dataset = dataset.map(parse_image, num_parallel_calls=5)

# Break your long tensors into many little ones
dataset = dataset.flat_map(
      lambda features: tf.data.Dataset.from_tensor_slices(features)
)

# Turn the dictionary in each record into a tuple without a label
dataset = dataset.map(
      lambda data_dict: (tf.transpose(list(data_dict.values())), )
)

# Turn each patch into a batch
dataset = dataset.batch(patch_width * patch_height)


# NOTE: assuming all data, regardless of file format or dataset,
#       have 5 bands (1 constant term (i.e. intercept),
#       and 1 per harmonic-term coefficient from the harmonic regressions)
bands = data.bands()

# get lists of lon, lat vals
lons = np.arange(data.ulx, data.lrx, data.resx)
lats = np.arange(data.uly, data.lry, data.resy)

# get 1 year of daily values, expressed in radians
annual_radian_days = np.linspace(0, 2*np.pi, 366)[:365]
semiannual_radian_days = np.linspace(0, 4*np.pi, 366)[:365] % (2 * np.pi)
sin1 = sin(annual_radian_days)
cos1 = cos(annual_radian_days)
sin2 = sin(semiannual_radian_days)
cos2 = cos(semiannual_radian_days)
harmonic_terms = [sin1, cos1, sin2, cos2]


#----------------
# run calculation
#----------------

# create the output asynch raster, to be filled in (starting as all NaNs)
asynch = beta_sin1*np.nan


# loop over pixels
for i, lat in enumerate(lats):
    for j, lon in enumerate(lons):

        # create lists of R2 and dist values
        R2s = []
        dists = []

        # calculate the focal pixel's time series
        ts_foc = calc_time_series(i, j, bands, harmonic_terms)

        # find the neighborhood of pixels for the focal pixel
        neigh_lats, neigh_lons = find_neighborhood_pixels(lat, lon)

        # loop over neighbors
        for ni, nlat in enumerate(neigh_lats):
            for nj, nlon in enumerate(neigh_lons):

                # get the neighbor's time series
                ts_neigh = calc_time_series(ni, nj, bands, harmonic_terms)

                # correlate the two time series and extract the R2
                R2 = linear_regression(ts_neigh, ts_foc)['R2']
                R2s.append(R2)

                # get the distance between the focal and neighbor pixels
                dist = calc_distance(lat, lon, nlat, nlon)
                dists.append(dist)
        
        # get the slope of the overall regression
        asynch_val = np.abs(linear_regression(R2s, dists)['slope'])
        output[i, j] = asynch_val


#------------
# data output
#------------

write_data(asynch, outfilepath)

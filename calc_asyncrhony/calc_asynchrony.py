#--------
# imports
#--------

import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import glob
import os
from pprint import pprint

#-----------
# set params
#-----------

MAX_DIST = 300000  # meters

DATA_DIR = ('/home/drew/Desktop/stuff/berk/research/projects/seasonality/'
              'GEE_output')

BANDS = ['constant', 'beta_sin1', 'beta_cos1', 'beta_sin2', 'beta_cos2']


#-----------------
# define functions
#-----------------

def calc_time_series(i, j, bands, harmonic_terms):
    """
    Calculates the time series at pixel i,j, using the coefficients for the
    sin and cosine terms of the annual and semiannual harmonic components.

    Returns the time series as a list.
    """
    constant, beta_sin1, beta_cos1, beta_sin2, beta_cos2 = bands
    sin1, cos1, sin2, cos2 = harmonic_terms
    ts = (constant[i, j] +
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
    pass


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


# data-parsing stuff from:
#        https://colab.research.google.com/github/google/
#        earthengine-api/blob/master/python/examples/ipynb/
#        TF_demo1_keras.ipynb#scrollTo=x2Q0g3fBj2kD)
def parse_tfrecord(example_proto):
    """
    The parsing function.
    
    Read a serialized example into the structure defined by featuresDict.

    Args:
        example_proto: a serialized Example.

    Returns:
        A tuple of the predictors dictionary and the label, cast to an `int32`.
    """
    parsed_features = tf.io.parse_single_example(example_proto, features_dict)
    pprint(parsed_features)
    labels = parsed_features.pop(LABEL)
    return parsed_features, tf.cast(labels, tf.int32)


def write_data(rast, filepath):
    """
    Writes the output raster to disk using the given filepath.
    """
    pass


#------------------
# handle input data
#------------------

# handle filepaths
infilepath = os.path.join(DATA_DIR, 'testRegCoeffExport-00000.tfrecord')
outfilepath = infilepath.split('.')[0]+'OUT.tfrecord'

# data-parsing stuff from:
#        https://colab.research.google.com/github/google/
#        earthengine-api/blob/master/python/examples/ipynb/
#        TF_demo1_keras.ipynb#scrollTo=x2Q0g3fBj2kD)
# Create a dataset from the TFRecord file in Cloud Storage.
dataset = tf.data.TFRecordDataset(infilepath)

# Get a list of fixed-length features, all of which are float32
columns = [
      tf.io.FixedLenFeature(shape=[1], dtype=tf.float32) for k in BANDS
]
features_dict = dict(zip(BANDS, columns))
pprint(features_dict)

# Map the data-parsing function over the dataset.
parsed_dataset = dataset.map(parse_tfrecord, num_parallel_calls=5)

# Print the first parsed record to check.
pprint(iter(parsed_dataset).next())

assert True == False

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

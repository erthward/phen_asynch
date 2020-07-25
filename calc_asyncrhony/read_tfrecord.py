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


#-----------
# set params
#-----------

DATA_DIR = ('/home/drew/Desktop/stuff/berk/research/projects/seasonality/'
              'GEE_output')
INFILEPATHS = glob.glob(os.path.join(DATA_DIR, '*.tfrecord'))

DIMS = (356, 356)

BANDS = ['beta_sin1', 'beta_cos1', 'beta_sin2', 'beta_cos2']


#--------------------
# get mixer file info
#--------------------

mixerfilepaths = glob.glob(os.path.join(DATA_DIR, '*mixer.json'))
assert len(mixerfilepaths) == 1, "MORE THAN 1 MIXER FILE FOUND!"
mixerfilepath = mixerfilepaths[0]

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


#-------------------
# read tfrecord file
#-------------------
raw_dataset = tf.data.TFRecordDataset(INFILEPATHS)

def parse_example_to_numpy(ex, dims, default_val=-9999.0):
    parsed = tf.train.Example.FromString(ex.numpy())
    feature = parsed.features.feature['constant']
    floats = feature.float_list.value
    arr = np.array(floats).reshape(dims)
    mask = np.int8(arr == default_val)
    print(mask)
    masked = np.ma.array(arr, mask=np.int8(arr == default_val))
    return masked

arrays = [parse_example_to_numpy(ex, DIMS) for ex in raw_dataset]
array_dict = dict(zip(BANDS, arrays))

print(array_dict)

fig = plt.figure()
n = 0
for k, v in array_dict.items():
    ax = fig.add_subplot(221+n)
    ax.imshow(v)
    ax.set_title(k)
    n+=1
plt.show()



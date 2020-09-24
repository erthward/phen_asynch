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

KERNEL_SIZE = 100

OVERLAP = int(KERNEL_SIZE / 2)

BANDS = ['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']


#-----------------
# define functions
#-----------------

def get_patch_dims(data_dir, kernel_size):
    """
    Gets the mixerfile located at data_dir.

    Returns the patch dimensions.
    """
    mixerfilepaths = glob.glob(os.path.join(data_dir, '*mixer.json'))
    assert len(mixerfilepaths) == 1, "MORE THAN 1 MIXER FILE FOUND!"
    mixerfilepath = mixerfilepaths[0]

    # read the mixer file
    mixer = json.load(open(mixerfilepath, 'r'))

    # Get relevant info from the JSON mixer file
    # NOTE: add overlap to account for the kernel size,
    # which isn't reflected in the mixer file
    patch_width = mixer['patchDimensions'][0] + kernel_size
    patch_height = mixer['patchDimensions'][1] + kernel_size
    patch_dimensions = (patch_width, patch_height)
    return(patch_dimensions)


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
        #mask = np.int8(arr == default_val)
        #print(mask)
        #masked = np.ma.array(arr, mask=np.int8(arr == default_val))
        arr[arr == default_val] = np.nan
        arrays.append(arr)
    out = np.stack(arrays)
    return out


def read_tfrecord_files(data_dir, kernel_size, bands):
    """
    Combines all the files contained in infilepaths into a TFRecordDataset.
    Returns a generator that yields each next example in the dataset, parsed
    into an n_bands x lon x lat numpy array.
    """
    # get the patch dimensions
    dims = get_patch_dims(data_dir, kernel_size)
    # get all the files to be read
    infilepaths = glob.glob(os.path.join(data_dir, '*.tfrecord'))
    # create a TFRecordDataset from the files
    raw_dataset = tf.data.TFRecordDataset(infilepaths)
    # turn the dataset into a list of arrays
    for example in raw_dataset.as_numpy_iterator():
        yield parse_tfexample_to_numpy(example, dims, bands)




#------------------
# TEST OUT THE CODE
#------------------

# get all 4 patches in the example files, concatenate into one map, and display
arrays = [arr for arr in read_tfrecord_files(DATA_DIR, KERNEL_SIZE, BANDS)]
concats = []

fig = plt.figure()
n = 0
for n, coeff in enumerate(BANDS):
    ax = fig.add_subplot(151+n)
    # concatenate the images, minus the overlap
    concat = np.concatenate([arr[n,
                                 OVERLAP:-OVERLAP,
                                 OVERLAP:-OVERLAP] for arr in arrays],
                           axis=0)
    med_val = np.nanmedian(concat)
    std_val = np.nanstd(concat)
    n_stds = 1
    ax.imshow(concat,
              vmin=med_val-n_stds*std_val,
              vmax=med_val+n_stds*std_val,
              cmap='twilight')
    ax.set_title(coeff)
    n+=1
    concats.append(concat)


# plot estimated SIF at new year's and mid-year for that 4-file concatenated
# dataset
fig2 = plt.figure()
ax1 = fig2.add_subplot(121)
ax2 = fig2.add_subplot(122)
jan_est = (concats[0] + concats[1]*np.sin(0) +
          concats[2]*np.cos(0) + concats[3]*np.sin(0) +
          concats[4]*np.cos(0))
ax1.imshow(jan_est, cmap='Spectral')
ax1.set_title('Jan. 1st')
jul_est = (concats[0] + concats[1]*np.sin(np.pi) +
          concats[2]*np.cos(np.pi) + concats[3]*np.sin(0) +
          concats[4]*np.cos(0))
ax2.imshow(jul_est, cmap='Spectral')
ax2.set_title('Jul. 1st')
plt.show()

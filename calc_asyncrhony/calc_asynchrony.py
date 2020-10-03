#!/bin/python
# calc_asynchrony.py

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

        In other words:

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

# TODO:
    # 0. figure out how to scale up premutation tests on GEE
    # 0.5 chat with Ian about asynch, R2s, ns, and scaling limitations
    # 1. figure out what to do about negative R2 values
    # 2. figure out how to write out the TFRecord files
    # 3. begin planning how to parallelize on XSEDE/Savio
    # 4. try to scale up the number of permutation tests on GEE
    # 5  figure out how to export the coeff patches globally, to Google Drive,
    #    then copy over to XSEDE/Savio from there

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

# directory where the data and mixerfile live
DATA_DIR = ('/home/drew/Desktop/stuff/berk/research/projects/seasonality/'
              'GEE_output')

# kernel size used by GEE to output the TFRecord files
KERNEL_SIZE = 100

# names of the bands saved into the TFRecord files
INBANDS = ['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']
OUTBANDS = ['asynch', 'asynch_R2', 'asynch_n']

# max distance out to which to find and include neighbors in
# each pixel's asynchrony calculation (in meters)
NEIGH_RAD = 150_000

# stdout options
VERBOSE = True
TIMEIT = True


#-----------------
# define functions
#-----------------

def get_infile_outfile_paths(data_dir):
    """
    Uses the data-directory path provided to get and return lists
    of the input and output files' paths.
    """
    # set up IO paths
    infilepaths = glob.glob(os.path.join(data_dir, '*.tfrecord'))
    outfilepaths = [fp.split('.')[0]+'OUT.tfrecord' for fp in infilepaths]
    return infilepaths, outfilepaths


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
    with tf.io.TFRecordWriter(filepath) as writer:
        for patch, label in zip(patches, labels):
            # serialize to bytes
            patch_bytes = tf.io.serialize_tensor(patch)
            # grab the shape
            patch_shape = patch.shape
            # serialize to a TF example
            example = serialize_tfrecord_example(patch_bytes, bands)
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
               patch[i, :, :].flatten())) for i, band in enumerate(bands)}
    #  create a Features message (i.e. protocol buffer) using tf.train.Example
    example_prot = tf.train.Example(features=tf.train.Features(feature=feature))
    return example_prot.SerializeToString()


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
    #gridx, gridy band [coords.flatten() for coords in np.meshgrid(xs, ys)]
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


def run_linear_regression(y, x, fit_intercept=True):
    """
    Calculates the simple linear regression of y on x.

    Returns the regression as a model object.
    """
    # Reshape the x variable
    reshaped_x = np.array(x).reshape(-1,1)
    # Make a model
    model = LinearRegression(fit_intercept=fit_intercept).fit(reshaped_x, y)
    # get the R2 and coeff values
    R2 = model.score(reshaped_x, y)
    slope = model.coef_ 
    intercept = model.intercept_
    return {'slope': slope, 'intercept': intercept, 'R2': R2}


def calc_distance(lon1, lat1, lon2, lat2):
    """
    Calculates the real-world distance (in meters)
    between the two points at lon1,lat1 and lon2,lat2.

    Returns the distance as a float.
    """

    pass


def get_neighbors_info(i, j, ys, xs, yres, xres, patch_dims,
                       tree=None, indices=None,
                       neigh_rad=NEIGH_RAD, min_dist_deg_lat=111_320,
                       deg_neigh_rad_buff_factor=1.05):
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
    (This actually gets way more pixels than necessary, because of all the
    pixels in the corners of the resulting box, which will be far outside
    the neighborhood radius. That said, it is a first-pass method for grabbing
    a set that contains all necessary pixels. 
    I SHOULD CONSIDER FINDING A MAX DECIMAL-DEGREE DISTANCE WITHIN WHICH TO
    GRAB PIXELS, THEN USING A kdTree INSTEAD!)

    Then it uses geopy's geodesic distance function to calculate the distance,
    in meters, of each pixel's centerpoint from the focal pixel's centerpoint.

    Only pixels whose centers are <= neigh_rad of the focal cell's center are
    retained, and those pixels' coordinates and distances are returned as keys
    and values in a dict.
    """
    if tree is None and indices is None:
        # get the minimum distance of one degree longitude at the maximum
        # absolute-value latitude in the input ys
        max_abs_lat = np.max(np.abs(ys))
        min_dist_deg_lon = geodesic((max_abs_lat, 0), (max_abs_lat, 1)).m

        # get number of x cells equiv to min_deg_lon_dist, to use as number of
        # neighbors surrounding the focal pixel in each direction that should be
        # included in the subsetted neighborhood raster
        # (cell_tot = m_tot / (m / deg) / (deg / cell)
        neigh_rad_y = int(np.ceil(np.abs(neigh_rad/min_dist_deg_lat/yres)))
        neigh_rad_x = int(np.ceil(np.abs(neigh_rad/min_dist_deg_lon/xres)))

        # use that number to get all the lats and lons from among which neighbors
        # will be chosen
        valid_ys = ys[i-neigh_rad_y:i+neigh_rad_y+1,
                      j-neigh_rad_x:j+neigh_rad_x+1].flatten()
        valid_xs = xs[i-neigh_rad_y:i+neigh_rad_y+1,
                      j-neigh_rad_x:j+neigh_rad_x+1].flatten()
        assert (valid_xs.size == valid_ys.size), ('Incorrect number of '
                #((neigh_rad_x*2)+1)*((neigh_rad_y*2)+1)), ('Incorrect number of '
                                                           'valid x and y coords '
                                                           ' returned!')
        # get the array indices corresponding to those coords
        valid_y_inds = np.indices(patch_dims)[0,:,:][i-neigh_rad_y:i+neigh_rad_y+1,
                      j-neigh_rad_x:j+neigh_rad_x+1].flatten()
        valid_x_inds = np.indices(patch_dims)[1,:,:][i-neigh_rad_y:i+neigh_rad_y+1,
                      j-neigh_rad_x:j+neigh_rad_x+1].flatten()
        assert (valid_y_inds.size == valid_x_inds.size == valid_xs.size), ('Number '
                   'of coord-indices not equal to number of coords!')
 
    else:
        # NOTE: buffer the max distance in degrees just a bit, so
        # that all points within the neighborhood distance in meters
        # are definitely found
        max_dist_deg = (neigh_rad/min_dist_deg_lat) * deg_neigh_rad_buff_factor
        # get the tree's data-row numbers for all neighbors
        # within the required distance
        neighs = tree.query_ball_point([xs[i, j], ys[i, j]], max_dist_deg)
        # use the row numbers to subset the pixels' centerpoint coords and
        # their index coordinate pairs
        valid_coords = tree.data[neighs, :]
        valid_xs = valid_coords[:,0]
        valid_ys = valid_coords[:,1]
        valid_indices = indices[neighs, :]
        valid_x_inds = valid_indices[:,0]
        valid_y_inds = valid_indices[:,1]
        assert (valid_xs.size == valid_ys.size), ('Incorrect number of '
            #((neigh_rad_x*2)+1)*((neigh_rad_y*2)+1)), ('Incorrect number of '
                                                       'valid x and y coords '
                                                       ' returned!')
        assert (valid_y_inds.size == valid_x_inds.size == valid_xs.size), ('Number '
               'of coord-indices not equal to number of coords!')
   
    # save the focal pixel's centerpoints, and lists of the other
    # centerpoints and their indices
    foc = (ys[i,j], xs[i,j])
    pts = [(valid_ys[n], valid_xs[n]) for n in range(len(valid_xs))]
    pts_inds = [(valid_y_inds[n], valid_x_inds[n]) for n in range(len(
                                                                valid_x_inds))]

    # calculate the geodesic distances
    # NOTE: per the Geopy docs
    # (https://geopy.readthedocs.io/en/latest/#module-geopy.distance)
    # I'm assuming the WGS84 ellipsoid model, and ignoring altitude differences
    # between points; seems perfectly fine as far as I'm concerned
    dists = [geodesic(foc, pt).m for pt in pts]
    out = {pt: [dist, pt_inds] for pt, dist, pt_inds in zip(pts, dists,
                                                           pts_inds) if (
                                                            dist <= neigh_rad)}
    return out


def make_cKDTree_and_indices(xs, ys):
    """
    Makes a cKDTree out of all the coordinates in a patch, to be used in
    neighbor-searching.
    """
    tree = cKDTree(np.stack((xs.flatten(), ys.flatten())).T)
    return tree


def make_indices_array(dims):
    """
    Makes an nx2 array of the n indices of the patch's pixels,
    where n = x_dim * y_dim.
    The structure of this indices array makes it subsettable by the output
    from a neighbor search on a cKDTree.
    """
    indices3D = np.indices(dims)
    x_indices = indices3D[1, :, :]
    y_indices = indices3D[0, :, :]
    indices = np.stack((x_indices.flatten(), y_indices.flatten())).T
    return indices


def calc_asynch_one_pixel(i, j, patch, patch_n, ys, xs, yres, xres, dims,
                          design_mat, tree, indices,
                          timeit=True, verbose=True):
    """
    Calculates and returns the asynchrony value for a single pixel
    located at position i,j in the patch.
    """
    if verbose and timeit:
        # get start time
        start = time.time()

    if verbose:
        print('\nPROCESSING PATCH %i: PIXEL (%i, %i)...\n' % (patch_n, i, j))
       
    # create lists of R2 and dist values
    R2s = []
    dists = []

    # calculate the focal pixel's time series
    ts_foc = calc_time_series(patch, i, j, design_mat)
 
    # get the coords, dists, and array-indices
    # of all of the focal pixel's neighbors
    coords_dists_inds = get_neighbors_info(i, j, ys, xs, yres,
                                           xres, dims,
                                           tree=tree,
                                           indices=indices)
    # loop over neighbors
    for neigh_coords, neigh_info in coords_dists_inds.items():
        # unpack the neighbor info
        neigh_dist, neigh_inds = neigh_info
        ni, nj = neigh_inds
    
        # get the neighbor's time series
        ts_neigh = calc_time_series(patch, ni, nj, design_mat)

        # drop this pixel if it returns NAs
        if np.any(np.isnan(ts_neigh)):
            pass
        else:
            # append the distance
            dists.append(neigh_dist)

            # calculate and append the R2
            R2 = run_linear_regression(ts_foc, ts_neigh)['R2']
            R2s.append(R2)
                   
    # get the slope of the overall regression
    # NOTE: setting fit_intercept to False and subtracting 1
    #       from the array of R2s effectively
    #       fixes the intercept at R2=1
    res = run_linear_regression(np.array(R2s) - 1,
                                    dists, fit_intercept=False)
    asynch = np.abs(res['slope'])
    asynch_R2 = res['R2']
    asynch_n = len(dists)

    if verbose and timeit:
        # get finish time
        stop = time.time()
        diff = stop-start
        print('\truntime: %0.4f' % diff)

    return asynch, asynch_R2, asynch_n


def calc_asynch(infilepath, outfilepath, inbands, outbands, dims,
                xmin, ymin, xres, yres, row_i, col_j, patch_n,
                verbose=True, timeit=True):
    """
    Read the patches from the input file, calculate and store
    the asynch metrics for each patch, then write the results
    to the output file.
    """
    # read the file's data in as a set of examples (i.e. patches)
    patches = read_tfrecord_file(infilepath, dims, inbands)

    # create a data container for all the asynch patches
    outpatches = []

    # loop over the patches within the current file
    # NOTE: each patch is of shape (n_bands, lat, lon)
    for patch in patches:

        # get the lons and lats of the current example's patch
        xs, ys = get_patch_lons_lats(xmin, ymin, xres, yres, dims, col_j, row_i)

        # get the cKDTree for these lons and lats
        tree = make_cKDTree_and_indices(xs, ys)

        # create the output asynch array, to be filled in (starts as all NaNs)
        # (and also an array to hold the R2s of the asynch-calc regressions)
        # (and also an array to hold the sample sizes for each of the asynch
        # regressions)
        asynch = np.nan * patch[0,:,:]
        asynch_R2s = np.nan * patch[0,:,:]
        asynch_ns = np.nan * patch[0, :, :]
        
        #----------------
        # run calculation
        #----------------

        # loop over pixels
        for i in range(patch.shape[1]):
            for j in range(patch.shape[2]):
 
                # leave the asynch output val as a NaN if the coeffs for this
                # pixel contain NaNs
                if np.any(np.isnan(patch[:, i, j])):
                    if verbose:
                        print(('\nPROCESSING PATCH %i: PIXEL (%i, %i)...'
                               '\n\n') % (patch_n, i, j))
                        print('\truntime: N/A')
                else:
                    (asynch_val,
                     asynch_R2_val,
                     asynch_n_val) = calc_asynch_one_pixel(i, j, patch, patch_n,
                                                           ys, xs, yres, xres,
                                                           dims, design_mat,
                                                           tree, indices,
                                                           verbose=verbose,
                                                           timeit=timeit)
                    asynch[i, j] = asynch_val
                    asynch_R2s[i, j] = asynch_R2_val
                    asynch_ns[i, j] = asynch_n_val
    
        # store this patch's asynch raster
        outpatch = np.stack((asynch, asynch_R2s, asynch_ns))
        outpatches.append(outpatch)

        #increment counters
        patch_n += 1
        if col_j == patches_per_row - 1:
            row_i += 1
            col_j = 0
        else:
            col_j += 1

        # TODO: DELETE ME
        #------------
        # data output
        #------------
        fig = plt.figure()
        ax1 = fig.add_subplot(221)
        ax1.set_title('asynch')
        ax1.imshow(asynch, cmap='magma')
        try:
            ax1.colorbar()
        except Exception:
            pass
        ax2 = fig.add_subplot(222)
        ax2.set_title('asynch_R2s')
        ax2.imshow(asynch_R2s, cmap='viridis')
        try:
            ax2.colorbar()
        except Exception:
            pass
        ax3 = fig.add_subplot(223)
        ax3.set_title('asynch_ns')
        ax3.imshow(asynch_ns, cmap='cividis')
        try:
            ax3.colorbar()
        except Exception:
            pass
        fig.show()
 
    
    # write out the asynch data
    write_tfrecord_file(outpatches, outfilepath, outbands)

#------------------------------------------------------------------------------
# calculate asynchrony
#---------------------

#---------------------
# set up IO paths, etc
#---------------------

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
# get the nx2 array of pixel index-pairs (which will be subsetted by cKDTree
# neighbor-query output
indices = make_indices_array(dims)

#--------------
# load the data
#--------------

# loop over the Tfrecord files, keeping track of the patch number and its
# corresponding row and column indices in the overall patch mosaic
patch_n = 0
row_i = 0
col_j = 0
for infilepath, outfilepath in zip(infilepaths, outfilepaths):
    calc_asynch(infilepath, outfilepath, INBANDS, OUTBANDS,
                dims, xmin, ymin, xres, yres, row_i, col_j, patch_n, 
                VERBOSE, TIMEIT)

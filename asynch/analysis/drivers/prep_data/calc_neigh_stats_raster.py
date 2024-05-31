#!/bin/python
# asynch_fns.py

"""
Reads in a global raster,
calculates and writes out mean, std, and/or shannon entropy
within pixel circular neighborhoods of a stipulated radius.
"""

#--------
# imports
#--------

from sklearn.linear_model import LinearRegression
from sklearn.neighbors import DistanceMetric
from sklearn.neighbors import BallTree
from haversine import haversine, Unit
from collections import Counter as C
from scipy.spatial import cKDTree
from scipy.stats import entropy
import matplotlib.pyplot as plt
from pprint import pprint
import geopandas as gpd
import rioxarray as rxr
import rasterio as rio
import xarray as xr
import numpy as np
import glob
import json
import time
import sys
import os


#-----------
# set params
#-----------

# directory where the data and mixerfile live
if os.path.abspath('.').split('/')[1] == 'home':
    data_dir = ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                'results/maps')
else:
    data_dir = '/global/scratch/users/drewhart/seasonality/'


# default missing-data val
default_vaL = -9999.0

# max distance out to which to find and include neighbors in
# each pixel's stat calculation
neigh_rad = 25_000

# convert that distance to degrees lat/lon, conservatively the ~47km/deg lon
# at 65 deg N/S that is the lowest value occurring within our dataset
neigh_rad_degs_conservative = neigh_rad/47_000

# stdout options
verbose = True
timeit = True

# get the input filename
filename = sys.argv[1]

# get stats to be calculated
assert len(sys.argv)>=3, 'Missing at least one stat to be calculated'
stats = sys.argv[2:]
for stat in stats:
    assert stat in ['m', 's', 'e'], ("Valid stats are 'm' (mean), "
                                     "'s' (std), "
                                     "and 'e' (Shannon entropy)")

# Shannon entropy function for vector of pixel values
def calc_shan_ent(vals):
    probs = [*{k: v/len(x) for k,v in C(x).items()}.values()]
    ent = entropy(probs)
    return ent

# dict of stats functions
stat_fns = {'e': calc_shan_ent,
            'm': np.mean,
            's': np.std,
           }

# dict of stats names
stat_names = {'e': 'Shannon entropy',
              'm': 'mean',
              's': 'standard deviation',
             }

#-----------------
# define functions
#-----------------


# function to find all neighbors of focal pixel within given radius
def calc_stats_one_pixel(x, y, da, out_da, neigh_rad,
                         neigh_rad_degs_conservative,
                         verbose=False, timeit=False):
    """
    1. use the conservative neighborhood radius (in degrees) to get the
    surrounding block of pixels that could be within the focal cell's
    neighborhood

    2. use haversine calculations to whittle that group down to the cells whose
    centroids actually are within the neighborhood

    3. return those cells' x,y coords
    """

    if (np.isnan(da.sel(x=x, y=y)) or
        np.allclose(da.attrs['_FillValue'], da.sel(x=x, y=y))):
        out_da.sel(x=x, y=y).values[0] = np.nan
    else:
        if verbose and timeit:
            start = time.time()

        # get the x and y resolutions of the raster
        xres = da.rio.transform().a
        yres = da.rio.transform().e

        # divide conservative degrees-based neighborhood radius by res, to get size
        # of potential neighborhood rectangle to subset
        x_rad_cells = neigh_rad_degs_conservative/xres
        y_rad_cells = neigh_rad_degs_conservative/yres

        # get the x and y slice edges
        x_min = x-(xres*x_rad_cells)
        x_max = x+(xres*x_rad_cells)
        # NOTE: y_min and y_max must be reversed because y_res is negative!
        assert yres < 1, "y_res not negative!"
        y_max = y-(yres*y_rad_cells)
        y_min = y+(yres*y_rad_cells)

        # grab the data subset in that rectangle
        # NOTE: this will automatically create less-than-full-size neighborhoods at
        #       raster edges, and then I loop over this sub_da's indices, so I
        #       don't need to manually check if x,y values are <=/>= the min x,y
        #       values in the dataset at any point
        sub_da = da.rio.clip_box(minx=x_min, miny=y_max, maxx=x_max, maxy=y_min)

        # find all cells actually within the neighborhood radius (in meters)
        # within that subset
        vals = []
        for nearby_x in sub_da.x:
            for nearby_y in sub_da.y:
                dist = haversine((y, x), (nearby_y, nearby_x), unit=Unit.METERS)
                if dist <= neigh_rad:
                    val = sub_da.sel(x=nearby_x, y=nearby_y).values[0]
                    vals.append(val)
        # append self's val to the list, too
        vals.append(sub_da.sel(x=x, y=y).values[0])

        # calculate and store statistics
        for i, stat in enumerate(stats):
            stat_val = stat_fns[stat](vals)
            if np.isinf(stat_val):
                stat_val = np.nan
            if verbose:
                print('\n\t%s = %0.4f\n\n' % (stat_names[stat], stat_val))
            out_da.sel(x=x, y=y).values[i] = stat_val
        if verbose and timeit:
            stop = time.time()
            timediff = stop-start
            print('\n\n\t\tpixel runtime: %0.1f seconds\n\n' % timediff)


def process_raster(filename, stats, neigh_rad, neigh_rad_degs_conservative,
                   verbose=False, timeit=False):
    # read the dataset
    da = rxr.open_rasterio(os.path.join(data_dir, filename))

    # set up coregistered output xr DataArray with one band for each stat to be
    # calculated
    out_da = xr.concat([da]*len(stats), dim=da.dims[0])
    # make all values NaN
    out_da[:,:,:] = np.nan

    # loop over x and y coordinates and fill up with neighborhood stats
    for x in da.x:
        for y in da.y:
            if verbose:
                print('\n\nNOW PROCESSING PIXEL AT (%0.4f, %0.4f)...\n\n' % (x, y))
            calc_stats_one_pixel(x, y, da, out_da, neigh_rad,
                                 neigh_rad_degs_conservative,
                                 verbose=verbose,
                                 timeit=timeit)

    # change DataArray attributes
    try:
        scale_factor = out_da.attrs['scale_factor']
    except Exception:
        scale_factor = 1.
    try:
        add_offset = out_da.attrs['add_offset']
    except Exception:
        add_offset = 0.
    try:
        _FillValue = out_da.attrs['_FillValue']
    except Exception:
        _FillValue = -9999
    out_da.attrs.clear()
    out_da.attrs['scale_factor'] = scale_factor
    out_da.attrs['add_offset'] = add_offset
    out_da.attrs['_FillValue'] = _FillValue
    out_da.attrs['long_name'] = stats

    # write to disk
    outfilename = os.path.join(data_dir, os.path.splitext(filename)[0] + 'STATS.tif')
    out_da.rio.to_raster(outfilename)


# process the requested data
if __name__ == "__main__":
    process_raster(filename, stats, neigh_rad, neigh_rad_degs_conservative,
                   verbose=verbose, timeit=timeit)

if verbose:
    print('\n%s\n\nstats calculationg for %s complete.' % ('-'*80, filename))


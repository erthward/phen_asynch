#!/usr/bin/python

# flake8: noqa

import os
import re
import datetime
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, MultiPolygon, Point


######
# TODO
######

    #1 decide whether I should be using 757nm or 771nm original SIF

    #2 instead of using all the polygons, just make a grid of lower-left grid
    #  corners, then subtract some amount from the centroid of each footprint
    #  and round to the closest cell in the grid, to find grid cells that are
    #  some effective buffer distance from any orbital track



###########
# FUNCTIONS
###########

def convert_iso_date(iso_date):
    date = datetime.date(*[int(x) for x in iso_date.split('-')])
    return(date)


def convert_netcdf_date(netcdf_date, filetype='orig'):
    if filetype == 'orig':
        date = '20%s-%s-%s' % (netcdf_date[:2], netcdf_date[2:4],
                               netcdf_date[4:6])
    elif filetype == 'grid':
        date = '%s-%s-%i' % (netcdf_date[:4],
                             netcdf_date[4:6],
                             (1 + (15 * netcdf_date[-1] == 'b')))
    elif filetype == 'trop':
        print('\n\nTODO: COMPLETE FOR TROPOMI')
        return
    else:
        print("\n'filetype' argument not valid!\n")
        return
    date = convert_iso_date(date)
    return(date)


def get_netcdf_files(data_dir='.', filetype='orig'):
    # create a regex pattern for SIF data filenames (which will capture the
    # year, month, and day from each file)
    if filetype == 'orig':
        patt = 'oco2_LtSIF_\d{6}_B8100r_\d*s\.nc4$'
        date_patt = patt = '(?<=oco2_LtSIF_)\d{6}(?=_B8100r_\d*s\.nc4)'
    elif filetype == 'grid':
        patt = 'sif_ann_\d{6}[ab]\.nc$'
        date_patt = patt = '(?<=sif_ann_)\d{6}[ab](?=\.nc)'
    elif filetype == 'trop':
        print('\n\nTODO: COMPLETE FOR TROPOMI')
        return
    else:
        print("\n'filetype' argument not valid!\n")
        return
    # create a list of the netCDFs
    files = [os.path.join(data_dir,
                          f) for f in os.listdir(
            data_dir) if re.search(patt, f) and not f.endswith('xml')]
    dates = [convert_netcdf_date(netcdf_date=re.search(
                        date_patt, f).group(),
                        filetype=filetype) for f in files]
    return(files, dates)


def read_netcdf(f, filetype='orig'):
    # check filetype arg
    if filetype not in ('orig', 'grid', 'trop'):
        print("\n'filetype' argument not valid!\n")
        return

    # get the appropriate null value
    null_val_dict = {'orig': -999999.0,
                     'grid': -999,
                    #TODO COMPLETE FOR TROPOMI
                     'trop': None
                    }
    null_val = null_val_dict[filetype]

    # get the appropriate variable names for lon, lat, and SIF vars
    netcdf_vars_dict = {'orig': ('footprint_vertex_longitude',
                                 'footprint_vertex_latitude',
                                 #'SIF_757nm'),
                                 'SIF_771nm'),
                        'grid': ('longitude',
                                 'latitude',
                                 'sif_ann'),
                        #TODO COMPLETE FOR TROPOMI
                        'trop': None
                       }


    # open dataset
    fh = Dataset(f, mode='r')
    # get the longitudes and  latitudes (which are 4-tuples
    # of the vertices of each footprint parallelogram)
    # and SIF values (for only 757nm in this case)
    lon_var, lat_var, sif_var = netcdf_vars_dict[filetype]
    lons = fh.variables[lon_var][:]
    lats = fh.variables[lat_var][:]
    sif = fh.variables[sif_var][:]

    # set the null value to that used in the netCDF

    # get lon_0 and lat_0 values to center the plot on
    #lon_0 = np.mean(plt_x_bnds)
    #lat_0 = np.mean(plt_y_bnds)
    # lat_0 = lats[np.where(lats != null_val)].mean()

    # close the netCDF file
    fh.close()

    return(lons, lats, sif)#, lon_0, lat_0)


def make_grid_cell_polygon(center_x, center_y, cell_size):
    diff = cell_size / 2
    xs = [center_x - diff, center_x + diff, center_x + diff, center_x - diff]
    ys = [center_y - diff] * 2 + [center_y + diff] * 2
    p = Polygon(zip(xs, ys))
    return p


# get the nearest value in an array either to the left ('down') or right ('up')
# of a certain value
def round_to_array(array, val, direction):
    if direction == 'down':
        comparison = idx = np.where(array < val)[0]
        if len(comparison) == 0:
            idx = 0
        else:
            idx = comparison.max()
    elif direction == 'up':
        comparison = idx = np.where(array > val)[0]
        if len(comparison) == 0:
            idx = array.size - 1
        else:
            idx = comparison.min()

    return idx


# get list of grid cells that have coverage within original OCO-2 orbital bands
#def get_data_coverage(orig_lons, orig_lats, grid_cells):
#    polys = [Polygon(zip(lon, lat)) for lon, lat in zip(orig_lons,
#                                                                  orig_lats)]
#    for p in polys:
#        if np.any(np.array([*p.exterior.coords.xy]) == -999999):
#            print('CONTAINS NULL:', p.exterior.coords.xy)
#        if not p.is_valid:
#            print('NOT VALID:', p.exterior.coords.xy)
#    mp = MultiPolygon([Polygon(zip(lon, lat)) for lon, lat in zip(orig_lons,
#                                                                  orig_lats)])
#    covered_cells = []
#    for n, cell in enumerate(grid_cells):
#        if mp.intersects(cell):
#            covered_cells.append(n)
#
#    return covered_cells



# get list of grid cells that have coverage within original OCO-2 orbital bands
def get_data_coverage(orig_lons, orig_lats,
                      grid_x_mins, grid_x_maxs,
                      grid_y_mins, grid_y_maxs,
                      coverage_array):
    #NOTE: max distance between parallelogram verices in either lat or lon
    #that I have observed is 0.11398697, so I'll buffer each parallelogram's
    #centroid by 0.15 degrees to be 'safe'
    buff_dist = 0.15
    # get an array of the parallelogram centroids
    centroids = np.vstack((np.mean(orig_lons, axis=1),
                           np.mean(orig_lats, axis=1))).T
    # for each centroid, get the max and min lon and lat indexes (in the
    # coverage array) that are overlapped by the coarse buffer added around the
    # pixel parallelogram, then use those indices to plop 1s into the coverage
    # array
    for cent_x, cent_y in centroids:
        print(cent_x, cent_y)
        min_x_idx = round_to_array(grid_x_mins, cent_x-buff_dist, 'down')
        max_x_idx = round_to_array(grid_x_maxs, cent_x+buff_dist, 'up')
        min_y_idx = round_to_array(grid_y_mins, cent_y-buff_dist, 'down')
        max_y_idx = round_to_array(grid_y_maxs, cent_y+buff_dist, 'up')
        print(min_x_idx, max_x_idx, min_y_idx, max_y_idx)
        coverage_array[min_y_idx:max_y_idx, min_x_idx:max_x_idx] = 1



###############################################
# FILES, VARIABLES, PARAMS, AND DATA STRUCTURES
###############################################

# get the data directory
data_dir = '/run/media/drew/SLAB/seasonality/SIF'
orig_dir = os.path.join(data_dir, 'OCO-2/orig')
grid_dir = os.path.join(data_dir,
                        'OCO-2/gridded/Global_High_Res_SIF_OCO2_1696/data')
trop_dir = os.path.join(data_dir, 'OCO-2/TROPOMI')

# set a date range to plot (written in ISO format)
start_date = '2016-01-01'
stop_date = '2017-12-31'

orig_files, orig_dates = get_netcdf_files(orig_dir)
grid_files, grid_dates = get_netcdf_files(grid_dir, filetype='grid')
# TODO GET TROPOMI FILES
#trop_files, trop_dates = get_netcdf_files(trop_dir)


# get the cell-centers from gridded ANN data
grid_lons, grid_lats, grid_sif_0 = read_netcdf(grid_files[0], filetype='grid')
grid_lon_mins = grid_lons - 0.025
grid_lon_maxs = grid_lons + 0.025
grid_lat_mins = grid_lats - 0.025
grid_lat_maxs = grid_lats + 0.025

#grid_lons_xi, grid_lats_yi = np.meshgrid(grid_lons, grid_lats)
#grid_cells = [make_grid_cell_polygon(lon,
#                                     lat, 0.05) for lon, lat in zip(
#                                grid_lons_xi.ravel(), grid_lats_yi.ravel())]

# create empty lists to fill up with all
#lon, lat coords for which there is original OCO-2 data
#covered_cells = set()

coverage_array = np.zeros((grid_lats.size, grid_lons.size))

#############################################################
# LOOP OVER AND READ IN FILES, COLLECTING GAP SITES FROM EACH
#############################################################

# now loop through all original OCO-2 files and determine which of those cells
# have no coverage, i.e. which cells are within the orbital gaps
for n, f in enumerate(orig_files):
    print('\nNow processing file number %i:\n\t%s\n\n' % (n, f))
    orig_lons, orig_lats, orig_sif = read_netcdf(f)
    print(len(orig_lons), len(orig_lats))
    # drop of footprints that seem to have missing vertex lat or lon vals
    # NOTE: not sure why this would be, but some parallelograms have
    # -999999 for some of their vertex values ...?!
    lons_null_rows = np.where(orig_lons == -999999)[0]
    lats_null_rows = np.where(orig_lons == -999999)[0]
    null_rows = [*set([*lons_null_rows] + [*lats_null_rows])]
    nonnull_rows = [row for row in range(
                            orig_lons.shape[0]) if row not in null_rows]
    orig_lons = orig_lons[nonnull_rows, :]
    orig_lats = orig_lats[nonnull_rows, :]
    print(len(orig_lons), len(orig_lats))

    get_data_coverage(orig_lons, orig_lats,
                      grid_lon_mins, grid_lon_maxs,
                      grid_lat_mins, grid_lat_maxs,
                      coverage_array)

    #out_cells = get_data_coverage(orig_lons, orig_lats, grid_cells)
    #covered_cells = covered_cells.union(set(out_cells))
    #print('Total number of covered cells:', len(covered_cells))
    print('Total number of covered cells:', np.sum(coverage_array == 1))
    print('---------------------')


# NEXT STEPS:

    #1. invert the coverage grid

    #2. use it and the grid_londs and grid_lats to select sample sites

    #3. extract both gridded and TROPOMI time series at those sample sites

    #4. figure out some way to reconcile times of those time series

    #5. run correlations and check strength

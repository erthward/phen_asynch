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

    #1 use square polygons instead of points for the grid cells

    #2 decide whether I should be using 757nm or 771nm original SIF



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


# get list of grid cells that have coverage within original OCO-2 orbital bands
def get_data_coverage(orig_lons, orig_lats, grid_lons, grid_lats):
    mp = MultiPolygon([Polygon(zip(lon, lat)) for lon, lat in zip(orig_lons,
                                                                  orig_lats)])
    covered_lons = []
    covered_lats = []
    for lon, lat in zip(grid_lons, grid_lats):
        if mp.contains(Point(lat, lon)):
            covered_lons.append(lon)
            covered_lats.append(lat)

    return covered_lons, covered_lats



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

# create empty lists to fill up with all lon, lat coords for which there is
# original OCO-2 data
covered_lons = convered_lats = set()

#############################################################
# LOOP OVER AND READ IN FILES, COLLECTING GAP SITES FROM EACH
#############################################################

# now loop through all original OCO-2 files and determine which of those cells
# have no coverage, i.e. which cells are within the orbital gaps
for n, f in enumerate(orig_files):
    if n < 2:
        print('\nNow processing file number %i:\n\t%s\n\n' % (n, f))
        orig_lons, orig_lats, orig_sif = read_netcdf(f)
        out_lons, out_lats = get_data_coverage(orig_lons, orig_lats,
                                               grid_lons, grid_lats)
        covered_lons = covered_lons.union(set(out_lons))
        covered_lats = covered_lats.union(set(out_lats))
        print(len(covered_lons))
        print(len(covered_lats))
        print('---------------------')











def rc():
    '''creates custom mapping settings, use mpl.rcdefaults() to default'''
    # can find out what is default by pprint(mpl.rcParams)
    mpl.rc('lines', linewidth=1.0, antialiased=True)  # didn't change
    mpl.rc('patch', linewidth=0.5, facecolor='348ABD', edgecolor='eeeeee',
           antialiased=True)
    mpl.rc('font', family='monospace', size=10.0)  # , weight="bold")
    mpl.rc('axes', facecolor='f4f4f4', edgecolor='black', linewidth=2,
           grid=False, titlesize='x-large', labelsize='x-large',
           labelcolor='black', axisbelow=True)
    # color_cycle='348ABD, 7A68A6, A60628, 467821, CF4457, 188487, E24A33')
    mpl.rc('xtick', color='black')
    mpl.rc('xtick.major', size=4, pad=6)
    mpl.rc('xtick.minor', size=2, pad=6)
    mpl.rc('ytick', color='black')
    mpl.rc('ytick.major', size=4, pad=6)
    mpl.rc('ytick.minor', size=2, pad=6)
    mpl.rc('legend', fancybox=True, fontsize=10.0)
    mpl.rc('figure', figsize='12, 7.5', dpi=88,
           facecolor='white')
    mpl.rc('figure.subplot', hspace=0.5,
           left=0.07, right=0.95, bottom=0.1,
           top=0.95)


def create_cmap():
    # make a red-to-green colormap
    colors = ['#980909', '#5C9809']
    cmap = LinearSegmentedColormap.from_list('my_map', colors, N=50)
    return(cmap)


def create_plot(rc, lat_0, lon_0):
    # create a figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # use custom matplotlib mappings defined above
    rc()

    # create a Basemap
    # m = Basemap(resolution='l',projection='cea',
    m = Basemap(resolution='l', projection='cyl',
                lat_ts=0, lat_0=lat_0, lon_0=lon_0, ax=ax)

    # create cmap
    cmap = create_cmap()

    # Add Grid Lines
    m.drawparallels(np.arange(-90, 90, 10.), labels=[1, 0, 0, 0], fontsize=8)
    merids = m.drawmeridians(np.arange(-180, 190, 10.), labels=[0, 0, 0, 1],
                             fontsize=8)

    # Add Coastlines, States, and Country Boundaries
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()

    # Add Title
    plt.title('SIF (OCO-2 data)')

    # Rotate the x labels
    for merid in merids.keys():
        try:
            merids[merid][1][0].set_rotation(45)
        except Exception:
            pass

    return(fig, ax, m, cmap)


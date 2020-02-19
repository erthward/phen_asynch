#!/usr/bin/python
# prelim_sif_seasonality_analysis.py

# flake8: noqa

########
# TODO:

########

import re
import os
import time
import datetime
import numpy as np
from netCDF4 import Dataset
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as mpl_Polygon
from matplotlib.colors import LinearSegmentedColormap


#################
# CHOOSE LOCATION
#################

# Bocas
plt_x_bnds = (-83.745346, -83.409576)
plt_y_bnds = (8.433659, 8.633297)

# Guajira Peninsula
plt_x_bnds = (-74.619141,-73.168945)
plt_y_bnds = (9.579084,10.336536)

# Manu National Park
plt_x_bnds = (-72.8167,-70.8175)
plt_y_bnds = (-12.1913,-10.6963)

# D.R. transect
#plt_x_bnds = (-71.065063,-70.101013)
#plt_y_bnds = (18.117140,20.012065)

# get D.R. transect 1km grid
#bbox = gpd.read_file('./DR_transect_1km_grid/DR_transect_1km_grid.shp')


def convert_netcdf_time(time):
    #NOTE: 17591 = March 01, 2018, according to the times provided
    #across the data files
    date = datetime.date(2018, 3, 1) + datetime.timedelta(days=int(
                                                                time) - 17591)
    return(date)


def get_netcdf_files(data_dir='.'):
    # create a regex pattern for SIF data filenames (which will capture the
    # year, month, and day from each file)
    patt = 'TROPO_SIF_\d{2}-2\d{3}\.nc$'
    date_patt = patt = '(?<=TROPO_SIF_)\d{2}-2\d{3}(?=\.nc)'
    # create a list of the netCDFs
    files = [os.path.join(data_dir,
                          f) for f in os.listdir(
            data_dir) if re.search(patt, f)]
    #dates = [convert_netcdf_time(time=re.search(
    #                    date_patt, f).group()) for f in files]
    #file_dict = dict(zip(dates, files))
    return(files)


def read_netcdf(f, null_val=-999.0):
    fh = Dataset(f, mode='r')

    # get the longitudes and  latitudes (which are 4-tuples
    # of the vertices of each footprint parallelogram)
    # and SIF values (for only 757nm in this case)
    lons = fh.variables['lon'][:]
    lats = fh.variables['lat'][:]
    times = fh.variables['time'][:]
    sif = fh.variables['dcSIF'][:]


    # subset for the chosen region
    keep_lons = np.bool8(
                np.int8(lons > plt_x_bnds[0]) * np.int8(lons < plt_x_bnds[1]))
    keep_lats = np.bool8(
                np.int8(lats > plt_y_bnds[0]) * np.int8(lats < plt_y_bnds[1]))
    keep_sif = sif[:, keep_lats, :][:, :, keep_lons]

    # get lon_0 and lat_0 values to center the plot on
    lon_0 = np.mean(plt_x_bnds)
    lat_0 = np.mean(plt_y_bnds)
    # lat_0 = lats[np.where(lats != null_val)].mean()

    # close the netCDF file
    fh.close()

    return(keep_lons, keep_lats, keep_sif, lon_0, lat_0, times)


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
    m.drawparallels(np.arange(-90, 90, 1), labels=[1, 0, 0, 0], fontsize=8)
    merids = m.drawmeridians(np.arange(-180, 190, 1), labels=[0, 0, 0, 1],
                             fontsize=8)

    # Add Coastlines, States, and Country Boundaries
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()

    # Add Title
    plt.title('SIF (TROPOMI data)')

    # Rotate the x labels
    for merid in merids.keys():
        try:
            merids[merid][1][0].set_rotation(45)
        except Exception:
            pass

    return(fig, ax, m, cmap)

def collect_regional_avg_data(sif, times):
    sif_means = np.mean(sif, axis = (1,2))
    data_dict = {convert_netcdf_time(times[i]): sif_means[i] for i in range(
                                                                len(times))}
    return data_dict


def plot_all_data_1cell(lat, lon):
    data_dict = {}
    # loop over monthly gridded files and pull out data
    for f in get_netcdf_files():
        fh = Dataset(f, mode='r')
        lons = fh.variables['lon'][:]
        lats = fh.variables['lat'][:]
        time = fh.variables['time'][:]
        sif = fh.variables['dcSIF'][:]
        time = [convert_netcdf_time(t) for t in time]
        # get the cell closest to the requested lat and lon 
        lat_idx = np.where(np.abs(lats - lat) == np.abs(lats - lat).min())
        lon_idx = np.where(np.abs(lons - lon) == np.abs(lons - lon).min())
        # add data to data_dict
        data_dict.update({time[n]: float(val) for n, val in enumerate(sif[:,
                                                                    lat_idx,
                                                                    lon_idx])})
    sorted_times = sorted([*data_dict.keys()])
    sorted_data = [data_dict[t] for t in sorted_times]
    plt.scatter(sorted_times, sorted_data, color='green')
    plt.show()
    return sorted_times, sorted_data

lat, lon = -8.553145, -64.705749
fig = plt.figure()
times, data = plot_all_data_1cell(lat, lon)



##################
# PROCESS THE DATA
##################

'''
# get all the files
files = get_netcdf_files()
print(files)
# files = [f for f in os.listdir() if f.startswith('TROPO_SIF')]

# create an empty data dictionary
data_dict = {}

# loop through files, reading in, averaging to the chosen area, and adding to
# the data dictionary
for f in files:
    lons, lats, sif, lon_0, lat_0, times = read_netcdf(f)
    data_dict.update(collect_regional_avg_data(sif, times))

# plot time series
xs = sorted([*data_dict.keys()])
ys = sorted([*data_dict.values()])
fig = plt.figure()
plt.suptitle('SIF (TROPOMI DATA)')
plt.xlabel('date')
plt.ylabel('SIF ($mW/m^2/sr/nm$)')
plt.plot(xs, ys)
plt.show()
'''

'''

# call all of the above functions to create my plot
files, dates = get_netcdf_files(data_dir)
start, stop = list(map(convert_iso_date, [start, stop]))
print('\nFOUND THE FOLLOWING FILES:\n%s' % '\n'.join(
    ['\t'+x for n, x in enumerate(files) if start <= dates[n] <= stop]))
print('----------------------------------------------------------\n\n\n')
plot_created = False
# save data for time-series plot
ts_dates = []
ts_data = []
gdfs = {}
for n, f in enumerate(files):
    # only plot data for the correct year (for now)
    if start <= dates[n] <= stop:
        print('\nPlotting date %s' % dates[n])
        lons, lats, SIF, lon_0, lat_0 = read_netcdf(f, null_val=-999999.0)
        if not plot_created:
            fig, ax, m, cmap = create_plot(rc=rc, lat_0=lat_0, lon_0=lon_0)
            plot_created = True
        ts_val, gdf = plot_data(ax=ax, cmap=cmap, lats=lats, lons=lons, data=SIF,
                           plot_map=plot_map)
        if not np.isnan(ts_val):
            ts_data.append(ts_val)
            ts_dates.append(dates[n])

            grid_overlay = gpd.overlay(bbox, gdf)
            gdfs[dates[n]] = grid_overlay


# show the plot
plt.show()

# plot the time-series data
ts_data_dict = dict(zip(ts_dates, ts_data))
fig2 = plt.figure()
plt.plot(sorted(ts_dates), [ts_data_dict[date] for date in sorted(ts_dates)])
plt.xlabel('date')
plt.ylabel('SIF value')
plt.show()

# TODO: choose a few cells across the DR transect (perhaps based on
# assessing how data rich each cell is), then plot the time series
# for those particular cells in a multi-line plot



##############################OLD CODE BELOW#########################
##############################OLD CODE BELOW#########################
##############################OLD CODE BELOW#########################
##############################OLD CODE BELOW#########################



def plot_data(ax, cmap, lats, lons, data, plot_map=True):
    # plot the footprint parallelograms that fall within
    # the requested range of longitudes (from plt_x_bnds and plt_y_bnds)
    # PLOT-COLOR APPROACH #1
    plot_data = []
    # PLOT-COLOR APPROACH #2
    patches = []
    colors = []

    # save data for time-series plot
    ts_data = []

    # create empty list to fill with geodataframes to concatenate
    gdfs = []

    for i in range(lons.shape[0]):
        f_p = Polygon(zip(lons[i, :], lats[i, :]))
        f_x, f_y = f_p.exterior.xy
        if ((np.int32(f_x) < plt_x_bnds[1]).all()
            and (np.int32(f_x) > plt_x_bnds[0]).all()
            and (np.int32(f_y) < plt_y_bnds[1]).all()
            and (np.int32(f_y) > plt_y_bnds[0]).all()):

            f_xi, f_yi = m(f_x, f_y)

            #save data for later cross-time analysis
            ts_data.append(SIF[i])
            geoseries = gpd.GeoSeries(f_p)
            series = pd.Series([SIF[i]])
            frame = {'sif':series, 'geo':geoseries}
            gdf = gpd.GeoDataFrame(frame, geometry='geo')
            gdfs.append(gdf)

            # PLOT-COLOR APPROACH #1
            polygon = mpl_Polygon(np.array(list(zip(f_xi, f_yi))), closed=True)
            patches.append(polygon)
            colors.append(cmap(SIF[i]))

            # PLOT-COLOR APPROACH #2
            plot_data.append(f_xi)
            plot_data.append(f_yi)
            plot_data.append(mpl.colors.rgb2hex(cmap(SIF[i])))
            # ax.fill(f_xi, f_yi, mpl.colors.rgb2hex(cmap(SIF[i])))

    # PLOT-COLOR APPROACH #1
    # ax.fill(*plot_data);

    # PLOT-COLOR APPROACH #2
    if plot_map:
        collection = PatchCollection(patches)
        ax.add_collection(collection)
        collection.set_color(colors)

    # average time-series data, concatenate gdfs, and return
    if len(gdfs) > 0:
        gdf = pd.concat(gdfs)
    else:
        gdf = None
    ts_data = np.mean(ts_data)
    return(ts_data, gdf)

'''

#!/usr/bin/python

# flake8: noqa

import os
import re
import sys
import datetime
import numpy as np
from netCDF4 import Dataset
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, MultiPolygon, Point
import statsmodels.api as sm

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf


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
                             (1 + (15 * (netcdf_date[-1] == 'b'))))
    elif filetype == 'trop':
        date = '%s-%s-01' % (netcdf_date[3:], netcdf_date[:2])
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
        patt = 'TROPO_SIF_\d{2}-2\d{3}\.nc$'
        date_patt = patt = '(?<=TROPO_SIF_)\d{2}-2\d{3}(?=\.nc)'
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

    # get the appropriate null/fill values (for lon, lat, and sif variables,
    # respectively)
    null_val_dict = {'orig': (9.969209968386869e+36,
                              9.969209968386869e+36,
                              #-999999.0),
                              9.969209968386869e+36),
                     'grid': (None,
                              None,
                              -999),
                     'trop': (9.969209968386869e+36,
                              9.969209968386869e+36,
                              -999.0)
                    }
    null_val = null_val_dict[filetype]

    # get the appropriate variable names for lon, lat, and SIF vars
    netcdf_vars_dict = {'orig': ('footprint_vertex_longitude',
                                 'footprint_vertex_latitude',
                                 'SIF_757nm'),
                                 #'SIF_771nm'),
                        'grid': ('longitude',
                                 'latitude',
                                 'sif_ann'),
                        'trop': ('lon',
                                 'lat',
                                 'dcSIF')
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
        #print(cent_x, cent_y)
        min_x_idx = round_to_array(grid_x_mins, cent_x-buff_dist, 'down')
        max_x_idx = round_to_array(grid_x_maxs, cent_x+buff_dist, 'up')
        min_y_idx = round_to_array(grid_y_mins, cent_y-buff_dist, 'down')
        max_y_idx = round_to_array(grid_y_maxs, cent_y+buff_dist, 'up')
        #print(min_x_idx, max_x_idx, min_y_idx, max_y_idx)
        coverage_array[min_y_idx:max_y_idx, min_x_idx:max_x_idx] = 1


def get_uncovered_sites_in_bbox(uncovered_lons, uncovered_lats,
                                llx, lly, urx, ury):
    out_lons = []
    out_lats = []
    for lon, lat in zip(uncovered_lons, uncovered_lats):
        if (llx <= lon <= urx) and (lly <= lat <= ury):
            out_lons.append(lon)
            out_lats.append(lat)

    return np.array((out_lons, out_lats)).T


def extract_vals_at_site(data, data_lons, data_lats, site_lon, site_lat,
                         filetype):
    lon = np.abs(data_lons - site_lon).argmin()
    lat = np.abs(data_lats - site_lat).argmin()
    if filetype == 'trop':
        extract = data[:, lat, lon]
    elif filetype == 'grid':
        extract = data[lat, lon]
    return extract


###############################################
# FILES, VARIABLES, PARAMS, AND DATA STRUCTURES
###############################################

# get the data directory
data_dir = '/media/deth/SLAB/diss/3-phn/SIF'
orig_dir = os.path.join(data_dir, 'OCO-2/orig')
grid_dir = os.path.join(data_dir,
                        'OCO-2/gridded/Global_High_Res_SIF_OCO2_1696/data')
trop_dir = os.path.join(data_dir,
                        'TROPOMI/fluo.gps.caltech.edu/data/tropomi/gridded')

# set a date range to plot (written in ISO format)
start_date = '2016-01-01'
stop_date = '2017-12-31'

orig_files, orig_dates = get_netcdf_files(orig_dir)
grid_files, grid_dates = get_netcdf_files(grid_dir, filetype='grid')
trop_files, trop_dates = get_netcdf_files(trop_dir, filetype='trop')

# get the cell-centers from gridded ANN data
grid_lons, grid_lats, grid_sif_0 = read_netcdf(grid_files[0], filetype='grid')
grid_lon_mins = grid_lons - 0.025
grid_lon_maxs = grid_lons + 0.025
grid_lat_mins = grid_lats - 0.025
grid_lat_maxs = grid_lats + 0.025

grid_lons_xi, grid_lats_yi = np.meshgrid(grid_lons, grid_lats)
#grid_cells = [make_grid_cell_polygon(lon,
#                                     lat, 0.05) for lon, lat in zip(
#                                grid_lons_xi.ravel(), grid_lats_yi.ravel())]

# create array of zeros, in which all cells with OCO-2
# coverage will be flipped to 1s
#coverage_array = np.zeros((grid_lats.size, grid_lons.size))
# NOTE: flipping because behavior of numpy or something changed at some point
#       and array starting reading in rightside up, but rest of code is written
#       for an array that was upside down!
coverage_array = np.flipud(np.loadtxt(('/media/deth/SLAB/diss/3-phn/SIF/OCO-2/'
                             'orbital_gap_validation_coverage_arrays/'
                             'final_coverage_array.txt')))

#############################################################
# LOOP OVER AND READ IN FILES, COLLECTING GAP SITES FROM EACH
#############################################################

# now loop through all original OCO-2 files and determine which of those cells
# have no coverage, i.e. which cells are within the orbital gaps
tmp_dir = '/media/deth/SLAB/diss/3-phn/SIF/OCO-2/orbital_gap_validation_coverage_arrays'

if 'final_coverage_array.txt' not in os.listdir(tmp_dir):
    # 02-08-2022: get rid of files already covered by existing tmp files
    # (array_up_to_file_XXXX.txt), so that I can rerun and pick up where
    # I left off if the previous job died midway
    last_file_completed = np.max([int(re.search(
                                    '(?<=^array_up_to_file_)\d+(?=.txt$)',
                                    f).group()) for f in os.listdir(tmp_dir)])
    orig_files_to_analyze = orig_files[last_file_completed+1:]

    for n, f in enumerate(orig_files_to_analyze):
        if n > -1:
            print('\nNow processing file number %i:\n\t%s\n\n' % (
                                                    n+last_file_completed+1, f))
            orig_lons, orig_lats, orig_sif = read_netcdf(f)
            #print(len(orig_lons), len(orig_lats))
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
            #print(len(orig_lons), len(orig_lats))

            get_data_coverage(orig_lons, orig_lats,
                              grid_lon_mins, grid_lon_maxs,
                              grid_lat_mins, grid_lat_maxs,
                              coverage_array)

            #out_cells = get_data_coverage(orig_lons, orig_lats, grid_cells)
            #covered_cells = covered_cells.union(set(out_cells))
            #print('Total number of covered cells:', len(covered_cells))
            print('Total number of covered cells:', np.sum(coverage_array == 1))
            print('---------------------')
            tmp_filename = 'array_up_to_file_%i.txt' % n
            if n % 5 == 0:
                np.savetxt(os.path.join(tmp_dir, tmp_filename), coverage_array,
                           fmt='%i')

else:
    print('\n\nCOVERAGE ARRAY ALREADY COMPLETE\n\n')


# get arrays of lon and lat values for the gridded data
grid_lon_i, grid_lat_j = np.meshgrid(grid_lons, grid_lats)

# get all lons and lats where the gridded data are not covered by orbitals
uncovered_lons = grid_lon_i[np.where(coverage_array == 0)]
uncovered_lats = grid_lat_j[np.where(coverage_array == 0)]

# get all uncovered lons and lats within a series of bounding boxes in tropical
# South America, Africa, and Australasia (Papua & Papua New Guinea, Java, NE
# Australia)
# NOTE: created all bounding boxes manually on
# https://boundingbox.klokantech.com
SA = -76.4851409407,-11.1511447176,-60.7356888937,7.6883728971
AF = 15.145367384,-26.9052129343,31.932476759,13.6216516376
PA = 138.7246979175,-6.3862012632,142.3686967951,-3.3381752281
JA = 108.1096467968,-7.6019672426,110.3972824892,-6.9865097164
AU = 141.8064223321,-19.8270873809,144.9807746735,-14.9180829021

# list to collect sampling sites for validation test
site_set = []
site_set_colors = []
site_set_color_dict = {SA: 'orange',
              AF: 'purple',
              PA: '#60c1cc', # 02-07-2022: changed all to same blue
              JA: '#60c1cc',
              AU: '#60c1cc',
             }
for n, box in [(60, SA), (60, AF), (20, PA), (20, JA), (20, AU)]:
    sites = get_uncovered_sites_in_bbox(uncovered_lons, uncovered_lats, *box)
    site_set.extend(sites[np.random.choice(range(sites.shape[0]), size=n,
                                           replace=False), :])
    site_set_colors.extend([site_set_color_dict[box]]*n)



# loop over TROPOMI and gridded SIF datasets, extracting sampling dates
# and values for each of our sites
trop_data_dates = []
trop_data_series = {tuple(site):[] for site in site_set}
for date, f in zip(trop_dates, trop_files):
    # save the data's date
    trop_data_dates.append(date)

    #open the file
    print('NOW EXTRACTING FROM: %s' % f)
    file_lons, file_lats, file_sif = read_netcdf(f, filetype='trop')
    #extract vals for each site and append to that site's list in the dict
    for site in site_set:
        # note: this gets all values for each day in the month
        site_data = extract_vals_at_site(file_sif, file_lons, file_lats,
                                         site[0], site[1], filetype='trop')
        trop_data_series[tuple(site)].append(site_data)


grid_data_dates = []
grid_data_series = {tuple(site):[] for site in site_set}
for date, f in zip(grid_dates, grid_files):
    if min(trop_data_dates) <= date <= datetime.date.fromordinal(
                                max(trop_data_dates).toordinal() + 30):
        #save the data's date
        grid_data_dates.append(date)

        #open the file
        print('NOW EXTRACTING FROM: %s' % f)
        file_lons, file_lats, file_sif = read_netcdf(f, filetype='grid')

        #extract val for each site and append it to that site's list
        for site in site_set:
            site_data = extract_vals_at_site(file_sif, file_lons, file_lats,
                                             site[0], site[1], filetype='grid')
            grid_data_series[tuple(site)].append(site_data)



# data for TROPOMI are a series of 28, 29, 30, or 31 values
# (or missing values) for a month, whereas gridded are a single value
# for either the 1st or 16th of a month
# so, go through the TROPOMI data and subset out just the values for the 1st
# and 16th of each month
trop_subset_data_series = {}
for site, series in trop_data_series.items():
    new_series = [*np.hstack([date_series[[0,15]] for date_series in series])]
    trop_subset_data_series[site] = new_series

# for each site, for each date, if the site has missing data in either dataset
# then set both dataset's values for that data to np.nan
for site, grid_series in grid_data_series.items():
    nans = np.where(np.array(grid_series) == -999)[0]
    masked = [i for i, n in enumerate(grid_series) if np.ma.is_masked(n)]
    missing = list(nans) + masked
    for idx in missing:
        grid_series[idx] = np.nan
        trop_subset_data_series[site][idx] = np.nan
for site, trop_series in trop_subset_data_series.items():
    missing = np.where(np.array(trop_series) == -999)[0]
    for idx in missing:
        trop_series[idx] = np.nan
        grid_data_series[site][idx] = np.nan


# run correlations and check their strength, potentially by region
corr_coeffs = {}
for site, grid_series in grid_data_series.items():
    trop_series = trop_subset_data_series[site]
    r = np.corrcoef([n for n in grid_series if not np.isnan(n)],
                    [n for n in trop_series if not np.isnan(n)])
    corr_coeffs[site] = r[0, 1]


####################
# PLOTS AND ANALYSIS
####################

# histogram of the correlation coefficients
fig_hist = plt.figure()
fig_hist.suptitle(('corr. coeffs.: TROPOMI vs. ANN-gridded OCO-2 '
                   'within orbital gaps'),
                  #fontdict={'fontsize': 25})
                 )
plt.hist([*corr_coeffs.values()], bins=25, alpha=0.8)
plt.xlabel('correlation coefficient')
plt.ylabel('frequency')


# line-plot grid for all sites' time series,
# and scatterplot of all TROPOMI-ANN pairs with regression results
def make_scatter_and_tsgrid():
    fig_scat = plt.figure(figsize=(20,10))
#    fig_scat.suptitle(('scatterplot of paired values across all sites, colored '
#                       'by region'),
#                      fontdict={'fontsize': 25})
    scat_ax = fig_scat.add_subplot(111)
    scat_ax.set_xlabel('ANN-gridded SIF ($mW/m^2/sr/nm$)',
                       fontdict={'fontsize':22})
    scat_ax.set_ylabel('TROPOMI SIF ($mW/m^2/sr/nm$)',
                       fontdict={'fontsize':22})
    # grid of paired time-series plots, one per site
    region_dict = {(-76.4851409407,-11.1511447176,
                    -60.7356888937,7.6883728971): 'S.Am.' ,
                   (15.145367384,-26.9052129343,
                    31.932476759,13.6216516376): 'Afr.',
                   (138.7246979175,-6.3862012632,
                    142.3686967951,-3.3381752281): 'Pap.',
                   (108.1096467968,-7.6019672426,
                    110.3972824892,-6.9865097164): 'Jav.',
                   (141.8064223321,-19.8270873809,
                    144.9807746735,-14.9180829021): 'Aus.'
                   }
    color_dict = {'S.Am.': 'orange',
                  'Afr.': 'purple',
                  'Pap.': '#60c1cc', # 02-07-2022: changed all to same blue
                  'Jav.': '#60c1cc',
                  'Aus.': '#60c1cc',
                 }
    fig_tsgrid = plt.figure()
    fig_tsgrid.suptitle(('paired TROPOMI (solid) and ANN-gridded (dashed) time '
                         'series, by site: ORANGE=South America, '
                         'PURPLE=Africa, BLUES=Papua New Guinea, Java, Australia '
                         '(light to dark)'))#, fontdict={'fontsize':25})
    fig_tsgrid.tight_layout()
    n_cols = 20
    n_rows = 9
    n = 1
    tot_series = [[], []]
    for site, grid_series in grid_data_series.items():
        #get the site's region
        region = None
        for reg, reg_name in region_dict.items():
            if ((reg[0] <= site[0] <= reg[2]) and
                (reg[1] <= site[1] <= reg[3])):
                region = reg_name
        # add an Axes object and plot the time series
        ax = fig_tsgrid.add_subplot(n_rows, n_cols, n)
        #ax.set_title(region + ': ' + str(site))
        ax.plot(grid_data_dates, trop_subset_data_series[site],
                '-', c=color_dict[region])
        ax.plot(grid_data_dates, grid_series, ':', c=color_dict[region])
        #increment my Axes counter
        n+=1
        # get rid of axis labels where necessary
        if n%20 != 2:
            ax.set_yticklabels([])
        else:
            plt.yticks(size=6)
            #ax.set_yticklabels(ax.get_yticklabels(), size=6)
        if n <= 161:
            ax.set_xticklabels([])
        else:
            plt.xticks(size=6, rotation=45)
            #ax.set_xticklabels(ax.get_xticklabels(), size=6, rotation=45)
        if n == 82:
            ax.set_ylabel('SIF ($mW/m^2/sr/nm$)', fontdict={'fontsize':22})
        if n == 171:
            ax.set_xlabel('time', fontdict={'fontsize':22})
        # add the site's TROPOMI and ANN-gridded data to the overall
        # scatterplot
        scat_ax.plot(grid_series, trop_subset_data_series[site], '.',
                      c=color_dict[region])
        # add the data to the tot_series lists
        tot_series[0].extend(grid_series)
        tot_series[1].extend(trop_subset_data_series[site])
    # calculate the linear regression between the two total series and plot it
    mod = sm.regression.linear_model.OLS(
                    endog = [n for n in tot_series[0] if not np.isnan(n)],
                    exog = [n for n in tot_series[1] if not np.isnan(n)]).fit()
    pred_xs = np.arange(0, 0.9, 0.01)
    pred_ys = mod.predict(pred_xs)
    scat_ax.plot(pred_ys, pred_xs, ':k')
    scat_ax.text(0.5, 0.25, '$R^2$: ' + str(np.round(mod.rsquared, 2)),
                 color='black', fontdict={'fontsize':20})
    scat_ax.text(0.5, 0.15, '$slope$: ' + str(np.round(mod.params[0], 2)),
                 color='black', fontdict={'fontsize':20})
    scat_ax.tick_params(labelsize=15)
    fig_scat.subplots_adjust(left=0.07, right=0.98, bottom=0.1, top=0.97)
    fig_scat.savefig(os.path.join(phf.FIGS_DIR,
                                  'FIG_SUPP_orbital_gap_validation_scatter.png'), dpi=500)
    fig_scat.show()
    fig_tsgrid.show()

# call that function
make_scatter_and_tsgrid()


def make_sample_point_map():
    fig_map = plt.figure(figsize=(20,10))
    map_ax = fig_map.add_subplot(111)
    cmap_list = [(0, '#ffffff'), (1, '#949c8e')]
    cmap = LinearSegmentedColormap.from_list('custom', cmap_list)
    map_ax.imshow(np.flipud(coverage_array), cmap=cmap, alpha=0.5)
    points = [(np.abs(grid_lons - site[0]).argmin(),
               np.abs(grid_lats[::-1]-site[1]).argmin()) for site in site_set]
    map_ax.scatter([p[0] for p in points], [p[1] for p in points],
                   c=site_set_colors, edgecolors='black', linewidth=0.4, s=8, alpha=0.7)
    map_ax.set_xticks([])
    map_ax.set_yticks([])
    fig_map.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)
    fig_map.show()
    fig_map.savefig(os.path.join(phf.FIGS_DIR,
                                 'FIG_SUPP_orbital_gap_validation_map.png'), dpi=500)


#call that function
make_sample_point_map()

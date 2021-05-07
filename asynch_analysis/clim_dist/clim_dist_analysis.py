#!/usr/bin/env python
# coding: utf-8

# py packages
import geopandas as gpd
import numpy as np
import glob
import json
import time
import os
#from osgeo import gdal
import random
import numpy as np
from scipy.spatial import KDTree
import math
import itertools
import rasterio as rio
import matplotlib.pyplot as plt
from shapely.ops import unary_union
from shapely.geometry import Polygon, Point
from scipy.spatial import ConvexHull
from collections import Counter as C

# local modules
import helper_fns as hf


# get the variables' filepaths
#BIO_DATA_DIR = ('C:\\Users\\thaon\\Documents\\Asynchrony\\bioclim')
DATA_DIR = ('/home/deth/Desktop/UCB/research/projects/seasonality/'
            'seasonal_asynchrony/asynch_analysis/data/')
BIOCLIM_DATA_DIR = DATA_DIR + 'bioclim/'
COUNTRIES_DATA_DIR = DATA_DIR + 'countries/'
COEFFS_FILE = DATA_DIR + 'SIF_coeffs.tif'
bioclim_infilepaths = glob.glob(os.path.join(BIOCLIM_DATA_DIR,"*.tif"))
bioclim_infilepaths = sorted(bioclim_infilepaths)

# load the country boundaries (to use for random point generation)
cntry = gpd.read_file(COUNTRIES_DATA_DIR + 'countries.shp')


#RANGEX = (-18000, 18000) # longtitude times 100 to get decimal value
#RANGEY = (-9000, 9000) # Latitude times 100 to get decimal value

n_pts = 1000  # Enter in the number greater than random points you need


# generate n random points within a given polygon
def generate_random_points_in_polygon(n, polygon):
    points = []
    minx, miny, maxx, maxy = polygon.bounds
    while len(points) < n:
        pnt = Point(random.uniform(minx, maxx),
                    random.uniform(miny, maxy))
        if polygon.contains(pnt):
            points.append(pnt)
    return points


##Generate random x,y coordinates
#def generate_random_coordinates(rangeX,rangeY, qty):
#    randPoints = [] # a list to store coordinates
#    while len(randPoints) < qty:
#        x = (random.randrange(*rangeX))/100 
#        y = (random.randrange(*rangeY))/100 
#        randPoints.append((x,y))
#    return randPoints
#
#
##Extract value at the coordinate
#def extract_bioclim_variables(randPoints):
#    list_bioclim = [] #a list to store bioclimatic variables
#    for point in randPoints:
#        list_bioclim += [point]
#        for path in bio_infilepaths:
#            raster = gdal.Open(path)
#            cols = raster.RasterXSize
#            rows = raster.RasterYSize
#            transform = raster.GetGeoTransform()
#            xOrigin = transform[0]
#            yOrigin = transform[3]
#            pixelWidth = transform[1]
#            pixelHeight = -transform[5]   
#            r = raster.GetRasterBand(1).ReadAsArray(0,0,cols,rows)
#            r = r.astype('float')
#            n = -3.4e+38 #value over the ocean to filter out
#            na = -3.3999999521443642e+38 #value over the ocean to filter out
#            r[r == z] = 'nan' 
#            r[r == na] = 'nan' 
#            r_mean = np.nanmean(r)
#            r_std = np.nanstd(r)
#            z_rast = (r - r_mean) / r_std
#            col = int((point[0] - xOrigin) / pixelWidth)
#            row = int((yOrigin - point[1] ) / pixelHeight)
#            list_bioclim +=[z_rast[row][col]]
#    return list_bioclim
#
#
## Functions return list of climate Euclidean distance between pairwise pixel
#def filter_out_ocean_pixel(biolist,coords, qty):
#    z = np.float32(-3.4e+38) #value of ocean region
#    bl = np.array_split(biolist, qty)
#    d = np.where(np.all(np.delete(bl,0,1) == np.repeat('nan',19),axis=1)) #index of array with pixel over the ocean
#    index = np.array(d).tolist() #index to delete
#    flat_list = [item for sublist in index for item in sublist]
#    for i in sorted(flat_list, reverse=True): #delete array with ocean pixel
#        del bl[i]
#    ncoords = coords
#    for j in sorted(flat_list, reverse=True):
#        del ncoords[j]
#    return bl, ncoords
#
#
#def calculate_euc_clim(biolist, coords):
#    pw = list(itertools.combinations(list(range(0,len(biolist))),2)) #generate list of pairwise combinations
#    ed = []
#    pw_coords = []
#    for i in pw:
#        dist = np.sqrt(np.sum((np.delete(biolist[i[0]],0) - np.delete(biolist[i[1]],0))**2))
#        ed += [dist]
#        pw_coords += [(coords[i[0]],coords[i[1]])]
#    return ed, pw_coords
#
#
#def calculate_euc_clim_at_random_points(rangeX,rangeY, qty): #Main Function
#    points = generate_random_coordinates(rangeX,rangeY, qty)
#    bioclim = extract_bioclim_variables(points)
#    bclim, bclim_coords = filter_out_ocean_pixel(bioclim,points, qty)
#    edbclim, bclim_pw_coords = calculate_euc_clim(bclim, bclim_coords)
#    return edbclim, bclim_pw_coords
#
#
#def get_seasonality_info_points(coeffs_rast_filepath, pts, dists=True):
#    """
#    takes a raster of coeffs from the harmonic seasonality regression
#    and an nx2 np.array of n points' x and y coordinates,
#    returns either the fitted time series at those points
#    (if 'dists' == False), or a matrix of pairwise seasonal distances
#    at those points (if 'dists' == True; default)
#    """
#    # read in the raster file
#    f = rio.open(coeffs_rast_filepath)
#    rast_coeffs = f.read()
#    #print(rast_coeffs.shape)
#
#    # get the cells' max lon- and lat-bound values
#    cell_max_lons, cell_max_lats = get_cell_lonlat_bounds(f)
#    #print('lons', cell_max_lons)
#    #print('lats', cell_max_lats)
#
#    # get the regression's design matrix
#    design_mat = make_design_matrix()
#
#    # matrix-multiply design_mat x coeffs to get fitted time series for all pts
#    ts_mat = np.zeros((pts.shape[0], 365)) * np.nan
#
#    for row_i in range(pts.shape[0]):
#        # get the array coords of the cell the point falls in
#        pt_cell_i, pt_cell_j = get_cell_coords_for_pt(lon=pts[row_i,0],
#                                                      lat=pts[row_i,1],
#                                                    cell_max_lons=cell_max_lons,
#                                                    cell_max_lats=cell_max_lats)
#        #print('pt_coeffs', rast_coeffs[:, pt_cell_i, pt_cell_j])
#        ts = np.sum(rast_coeffs[:, pt_cell_i, pt_cell_j] * design_mat, axis=1)
#        ts_mat[row_i, :] = ts
#
#    if dists:
#        # calculate pairwise Euc dist matrix
#        pw_dist_mat = np.zeros((pts.shape[0], pts.shape[0])) * np.nan
#        for row_i in range(pw_dist_mat.shape[0]):
#            for col_j in range(pw_dist_mat.shape[1]):
#                if row_i == col_j:
#                    dist = 0
#                else:
#                    dist = calc_euc_dist(ts_mat[row_i,:], ts_mat[col_j, :])
#                pw_dist_mat[row_i, col_j] = dist
#        return pw_dist_mat
#    else:
#        # return the time series matrix, if pairwise distances not requested
#        return ts_mat
#
#
#def get_cell_lonlat_bounds(rio_file):
#    """
#    takes a rasterio opened-file object,
#    returns a tuple of two arrays, the first for lons, the second for lats,
#    each of which contains a cell's max coordinate values for each col (lons)
#    or row (lats) in the raster
#    """
#    # get the res, min, and length values for lon and lat dimensions
#    aff = rio_file.transform
#    x_res = aff.a
#    x_min = aff.c
#    x_len = rio_file.width
#    y_res = aff.e
#    y_min = aff.f
#    y_len = rio_file.height
#    # get arrays of the max (i.e. east and south) lon and lat cell boundaries
#    # for all the cells in the raster
#    lons_east_cell_bounds = np.linspace(start=x_min+x_res,
#                                        stop=(x_res*x_len)+(x_min+x_res),
#                                        num=x_len)
#    lats_south_cell_bounds = np.linspace(start=y_min+y_res,
#                                         stop=(y_res*y_len)+(y_min+y_res),
#                                         num=y_len)
#    return lons_east_cell_bounds, lats_south_cell_bounds
#
#
#def get_cell_coords_for_pt(lon, lat, cell_max_lons, cell_max_lats):
#    """
#    uses provided dictionary to crosswalk a point's coordinates
#    with its cell's row,col (i.e. i,j) coordinates in the numpy array
#    holding a raster's data
#    """
#    #londiff_gtz = (cell_max_lons - lon) > 0
#    #latdiff_gtz = (cell_max_lats - lat) > 0
#    #j = np.argmax(~londiff_gtz[:-1] == londiff_gtz[1:]) + 1
#    #i = np.argmax(~latdiff_gtz[:-1] == latdiff_gtz[1:]) + 1
#    j = np.where((cell_max_lons - lon) > 0)[0][0]
#    i = np.where((cell_max_lats - lat) < 0)[0][0]
#    return i,j
#
#
#def calc_euc_dist(a1, a2):
#    """
#    Calculates the Euclidean distance between two 1d, length-n numpy arrays.
#    Returns the distance as a float.
#    """
#    dist = np.sqrt(np.sum((a1 - a2)**2))
#    return dist
#
#
#
##Sample code
#coeffs_rast_filepath = 'C:\\Users\\thaon\\Documents\\Asynchrony\\SIF_coeffs.tif'
#edist = []
#for i in range(0,len(tropicsample_coords)):
#    m = get_seasonality_info_points(coeffs_rast_filepath, np.array(pw_coords[i]), dists=True) #how should i store my points?
#    edist += [m[0][1]]
#
#
##Main Function to generate a graph. But for some reason, it cannot generate a graph more than 50 sample points. 
#def generate_graph_of_climate_asychrony(rangeX, rangeY, qty,coeffs_rast_filepath):
#    regionsample_bio, regionsample_coords = calculate_euc_clim_at_random_points(rangeX, rangeY, qty)
#    edist = []
#    for i in range(0,len(regionsample_coords)):
#        m = get_seasonality_info_points(coeffs_rast_filepath, np.array(regionsample_coords[i]), dists=True) #how should i store my points?
#        edist += [m[0][1]]
#    plt.scatter(regionsample_bio,edist)
#    plt.xlabel("Climate Euclidean Distance")
#    plt.ylabel("SIF Euclidean Distance")
#
#
##Tropical Region
#Tropic_RANGEX = (-18000, 18000) # longtitude
#Tropic_RANGEY = (-2000, 2000) # Latitude
#QTY = 30  # Enter in the number greater than random points you need
##tropicsample_bio, tropicsample_coords = calculate_euc_clim_at_random_points(Tropic_RANGEX,Tropic_RANGEY, QTY)
#generate_graph_of_climate_asychrony(Tropic_RANGEX,Tropic_RANGEY, QTY ,coeffs_rast_filepath)
#
#
##Sample plot
#plt.scatter(tropicsample_bio,edist)
#plt.xlabel("Climate Euclidean Distance")
#plt.ylabel("SIF Euclidean Distance")


def get_bioclim_dists_pts(pts, nodata_val=-3.4e+38):
    # list to hold arrays of all bioclim vars' vals
    bioclim_vals_arrs = []
    for fn in bioclim_infilepaths:
        # extract vals for each file
        bioclim_vals = hf.get_raster_info_points(fn, pts, 'vals')
        # mask out missing vals
        bioclim_vals[np.isclose(bioclim_vals, nodata_val)] = np.nan
        # standardize the vals
        stand_vals = (bioclim_vals -
                      np.nanmean(bioclim_vals))/np.nanstd(bioclim_vals)
        bioclim_vals_arrs.append(stand_vals)
    # concatenate into one array
    bioclim = np.concatenate(bioclim_vals_arrs, axis=1)
    # calculate pairwise dists
    pwd = np.zeros([pts.shape[0]]*2) * np.nan
    for i in range(pts.shape[0]):
        for j in range(i, pts.shape[0]):
            if i == j:
                pwd[i,j] = 0
            else:
                dist = hf.calc_euc_dist(bioclim[i,:],
                                        bioclim[j,:])
                pwd[i,j] = dist
                pwd[j,i] = dist
    return pwd



# define regions for:
# high asynchrony:
    # temperate: CA and SW
ha_tmp = Polygon([[-123.588671875,42.1726562502269],
                  [-120.951953125,35.75664368702668],
                  [-109.43828125,31.96275411512188],
                  [-105.746875,32.03728956520954],
                  [-106.362109375,43.653059057103164],
                  [-123.588671875,42.1726562502269]
                 ])
    # tropical: N Andes and N SA
ha_trp = Polygon([[-76.84073042457382,-0.40960685149885684],
                  [-58.12002729957382,-0.40960685149885684],
                  [-58.12002729957382,11.466203915130043],
                  [-76.84073042457382,11.466203915130043],
                  [-76.84073042457382,-0.40960685149885684]
                 ])
#ha_trp = Polygon([[-76,8.5],
#                  [-76,-1],
#                  [-58,-1],
#                  [-58,8.5],
#                  [-76,8.5]
#                 ])
# low asynchrony:
    # temperate: US South and SE
la_tmp = Polygon([[-95.39198282257864,30.76200580067661],
                  [-81.68104532257864,30.76200580067661],
                  [-81.68104532257864,37.797946552183234],
                  [-95.39198282257864,37.797946552183234],
                  [-95.39198282257864,30.76200580067661]
                 ])
    # tropical: W Amazon
la_trp = Polygon([[-72.70987104957382,-11.49384470495583],
                  [-58.12002729957382,-11.49384470495583],
                  [-58.12002729957382,-1.9912406261155682],
                  [-72.70987104957382,-1.9912406261155682],
                  [-72.70987104957382,-11.49384470495583]
                 ])
regs = [ha_tmp, ha_trp, la_tmp, la_trp]

# draw random points for each region
regs_pts = [generate_random_points_in_polygon(n_pts, reg) for reg in regs]

# create the figure
fig, ax = plt.subplots(1,1)
fig.suptitle(('seasonal distance vs. climatic distance, '
              'across latitude and asynchrony'), size=20)

# list of region names
reg_names = ['high asynch: temperate',
             'high asynch: tropical',
             'low asynch: temperate',
             'low asynch: tropical'
            ]

# create colors for the 4 regions
reg_cols = ["#7428ed", # ha_tmp
            "#18f04a", # ha_trp
            "#a6a4b3", # la_tmp
            "#8db08f", # la_trp
           ]

# dict to store the distances
dist_dict = {}

# get the nodata val
nodata_val = rio.open(bioclim_infilepaths[0]).nodata

for reg_poly, reg_pts, reg_col, reg_name in zip(regs,
                                                regs_pts,
                                                reg_cols,
                                                reg_names):

    print('\ngetting distance matrices for region: %s\n' % reg_name)

    # get points as nx2 numpy array
    pts = np.concatenate([np.array(pt.coords) for pt in reg_pts], axis=0)

    # get points' pairwise clim dists
    clim_dist = get_bioclim_dists_pts(pts, nodata_val=nodata_val)

    # get points' pairwise ts dists
    seas_dist = hf.get_raster_info_points(COEFFS_FILE, pts, 'ts_pdm')

    # drop clim dists for points without ts dists, and vice versa
    not_missing = np.where(np.nansum(seas_dist, axis=0)>0)[0]
    seas_dist = seas_dist[:, not_missing][not_missing,:]
    clim_dist = clim_dist[:, not_missing][not_missing,:]
    still_not_missing = np.where(np.nansum(clim_dist, axis=0)>0)[0]
    seas_dist = seas_dist[:, still_not_missing][still_not_missing,:]
    clim_dist = clim_dist[:, still_not_missing][still_not_missing,:]

    print(('\n%i points remain after dropping locations without seasonality '
          'data\n' % seas_dist.shape[0]))

    # extract the lower triangular values and scatter them
    indices = np.tril_indices(seas_dist.shape[0])

    seas_dist_vals = seas_dist[indices]
    clim_dist_vals = clim_dist[indices]

    # scatter 
    ax.scatter(clim_dist_vals, seas_dist_vals, s=1,
               c=reg_col, alpha=0.2, label=reg_name)

    # add the convex hull
    hull_pts = np.array((clim_dist_vals, seas_dist_vals)).T
    hull = ConvexHull(hull_pts)
    for simplex in hull.simplices:
        ax.plot(hull_pts[simplex, 0], hull_pts[simplex, 1], color=reg_col)

    # store the dists
    dist_dict[reg_name] = {'clim': clim_dist_vals,
                           'seas': seas_dist_vals
                          }

ax.legend(fontsize=16)

ax.set_xlabel('climate distance (Euclidean)', size=18)
ax.set_ylabel('seasonal distance (Euclidean)', size=18)
ax.tick_params(labelsize=16)
fig.show()

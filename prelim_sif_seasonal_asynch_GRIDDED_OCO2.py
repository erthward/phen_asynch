import re
import os
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


data_dir = ('/run/media/drew/SLAB/SIF_data/OCO-2/gridded/'
            'Global_High_Res_SIF_OCO2_1696/data')

files = {f: os.path.join(data_dir, f) for f in os.listdir(data_dir)}


# stipulate indices of lons and lats arrays corresponding to a few test sites
CA_coords = (1218, 2552)  # Californian Sierra Nevada: W 119.07499999998615, N 37.62499999999275
MX_coords = (1674, 2150)  # Mexican Sierra Madre Oriental: W 96.27499999998096, N 17.524999999993895
EC_coords = (2059, 1779)  # Ecuadorian Andes:  W 78.12499999997684, S 0.8750000000050591
#EC_coords = (2037, 1782)  # Ecuadorian Andes:  W 78.12499999997684, S 0.8750000000050591

def convert_iso_date(iso_date):
    date = datetime.date(*[int(x) for x in iso_date.split('-')])
    return date


def convert_netcdf_date(netcdf_date):
    date = '20%s-%s-%s' % (netcdf_date[:2], netcdf_date[2:4], netcdf_date[4:6])
    date = convert_iso_date(date)
    return date


def calc_dist(i, j, cent_ind):
    dx = np.abs(j - cent_ind)
    dy = np.abs(i - cent_ind)
    dist = np.sqrt(dx**2 + dy**2)
    return dist


def calc_R(x, y):
    corr_mat = np.corrcoef(x, y)
    R = corr_mat[0,1]
    return R


def read_netcdf(f, null_val=-999, max_neigh_dist=8):
    fh = Dataset(f, mode='r')

    # get the longitudes and  latitudes (which are 4-tuples
    # of the vertices of each footprint parallelogram)
    # and SIF values (for only 757nm in this case)
    lons = fh.variables['longitude'][:]
    lats = fh.variables['latitude'][:]
    # SIF = fh.variables['SIF_757nm'][:]
    sif = fh.variables['sif_ann'][:]

    # set the null value to that used in the netCDF

    # get lon_0 and lat_0 values to center the plot on
    lat_0 = lats[np.where(lats != null_val)].mean()
    lon_0 = lons[np.where(lons != null_val)].mean()

    # close the netCDF file
    fh.close()

    # get SIF vals for each site
    CA.append(sif[CA_coords[1], CA_coords[0]])
    MX.append(sif[MX_coords[1], MX_coords[0]])
    EC.append(sif[EC_coords[1], EC_coords[0]])

    #append rasts to raster-lists, for asynch calculation
    CA_rast = sif[CA_coords[1] - max_neigh_dist: CA_coords[1] + max_neigh_dist + 1,
                  CA_coords[0] - max_neigh_dist: CA_coords[0] + max_neigh_dist + 1]
    CA_asynch.append(CA_rast)
    MX_rast = sif[MX_coords[1] - max_neigh_dist: MX_coords[1] + max_neigh_dist + 1,
                  MX_coords[0] - max_neigh_dist: MX_coords[0] + max_neigh_dist + 1]
    MX_asynch.append(MX_rast)
    EC_rast = sif[EC_coords[1] - max_neigh_dist: EC_coords[1] + max_neigh_dist + 1,
                  EC_coords[0] - max_neigh_dist: EC_coords[0] + max_neigh_dist + 1]
    EC_asynch.append(EC_rast)

    # get date of dataset
    bn = os.path.split(f)[1]
    yr = bn[8:12]
    mn = bn[12:14]
    if bn[14] == 'a':
        dy = '1'
    elif bn[14] == 'b':
        dy = '16'
    else:
        raise ValueError

    iso_date = '%s-%s-%s' % (yr, mn, dy)
    date = convert_iso_date(iso_date)
    dates.append(date)

    return


# lists to hold focal site's sif vals
CA = []
MX = []
EC = []
dates = []


# lists to hold to neighborhood rasts, for asynch calculation
max_neigh_dist = 5
CA_asynch = []
MX_asynch = []
EC_asynch = []


for f in [*files.values()][:16]:
    read_netcdf(f, max_neigh_dist=max_neigh_dist)


# plot time series for the three sites
df = pd.DataFrame.from_dict({'date': dates,
                             'CA': CA,
                             'MX': MX,
                             'EC': EC})
    
df = df.sort_values(by='date')

fig = plt.figure()
plt.plot(df.date, df.CA, label = 'CA')
plt.plot(df.date, df.MX, label='MX')
plt.plot(df.date, df.EC, label='EC')
plt.legend()
plt.show()


#################################
# calc asynch for the three sites
#################################

# sort the raster-lists
sort_inds = []
date_dict = dict(zip(dates, range(len(dates))))
for d in sorted(dates):
    for k, v in date_dict.items():
        if k == d:
            sort_inds.append(v)

CA_asynch_sort = [CA_asynch[i] for i in sort_inds]
MX_asynch_sort = [MX_asynch[i] for i in sort_inds]
EC_asynch_sort = [EC_asynch[i] for i in sort_inds]

CA_stack = np.ma.dstack(CA_asynch_sort)
MX_stack = np.ma.dstack(MX_asynch_sort)
EC_stack = np.ma.dstack(EC_asynch_sort)

stack_dict = {'CA': CA_stack,
              'MX': MX_stack,
              'EC': EC_stack}


fig2 = plt.figure()
markers = ['o', 'x', 's']
Rs_dict = {}
dists_dict = {}
for n, item in enumerate(stack_dict.items()):
    loc, stack = item
    Rs = []
    dists = []
    cent_ind = int((stack.shape[0] - 1) / 2)
    focal_ts = stack[cent_ind, cent_ind, :]
    for i in [n for n in range(stack.shape[0]) if not n == cent_ind]:
        for j in [n for n in range(stack.shape[1]) if not n == cent_ind]:
            #TODO: calc distance in real distance units
            dists.append(calc_dist(i, j, cent_ind))
            Rs.append(calc_R(focal_ts, stack[i, j, :]))
    dists_dict[loc] = dists
    Rs_dict[loc] = Rs
    plt.scatter(dists, Rs, label=loc, alpha=0.5, marker=markers[n])
plt.legend()
plt.xlabel('neighbor distance')
plt.ylabel('R')
plt.show()

CA_res = np.corrcoef([n for n in Rs_dict['CA'] if not np.isnan(
            n)], [n for i, n in enumerate(dists_dict['CA']) if not np.isnan(
            Rs_dict['CA'][i])])[0,1] ** 2

MX_res = np.corrcoef([n for n in Rs_dict['MX'] if not np.isnan(
            n)], [n for i, n in enumerate(dists_dict['MX']) if not np.isnan(
            Rs_dict['MX'][i])])[0,1] ** 2

EC_res = np.corrcoef([n for n in Rs_dict['EC'] if not np.isnan(
            n)], [n for i, n in enumerate(dists_dict['EC']) if not np.isnan(
            Rs_dict['EC'][i])])[0,1] ** 2

print('CA result: %0.3f' % CA_res)
print('MX result: %0.3f' % MX_res)
print('EC result: %0.3f' % EC_res)

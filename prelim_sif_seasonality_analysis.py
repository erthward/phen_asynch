#!/usr/bin/python
#prelim_sif_seasonality_analysis.py

import re, os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as mpl_Polygon
from matplotlib.colors import LinearSegmentedColormap

#get the data directory
data_dir = '/media/ihavehands/SLAB/SIF_data/'

#set a range of longitudes within which to plot any data
plt_x = (-180, 180)

#set a date range to plot (written in ISO format)
start = '2015-05-07'
stop = '2015-05-08'


def convert_iso_date(iso_date):
    date = datetime.date(*[int(x) for x in iso_date.split('-')]) 
    return(date)


def convert_netcdf_date(netcdf_date):
    date = '20%s-%s-%s' % (netcdf_date[:2], netcdf_date[2:4], netcdf_date[4:6])
    date = convert_iso_date(date)
    return(date)


def get_netcdf_files(data_dir='.'):
    #create a regex pattern for SIF data filenames (which will capture the year, month, and day from each file)
    patt = 'oco2_LtSIF_\d{6}_B8100r_\d*s\.nc4$'
    date_patt = patt = '(?<=oco2_LtSIF_)\d{6}(?=_B8100r_\d*s\.nc4)'
    #create a list of the netCDFs
    files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if re.search(patt, f) and not f.endswith('xml')]
    dates = [convert_netcdf_date(netcdf_date = re.search(date_patt, f).group()) for f in files]
    return(files, dates)


def read_netcdf(f, null_val = -999999.0):
    fh = Dataset(f, mode = 'r')

    #get the longitudes and  latitudes (which are 4-tuples of the vertices of each footprint parallelogram)
    #and SIF values (for only 757nm in this case)
    lons = fh.variables['footprint_vertex_longitude'][:]
    lats = fh.variables['footprint_vertex_latitude'][:]
    #SIF = fh.variables['SIF_757nm'][:]
    SIF = fh.variables['SIF_771nm'][:]
    
    #set the null value to that used in the netCDF
    
    #get lon_0 and lat_0 values to center the plot on
    lon_0 = mean(plt_x)
    lat_0 = lats[np.where(lats != null_val)].mean()

    #close the netCDF file
    fh.close()

    return(lons, lats, SIF, lon_0, lat_0)


def rc():
    '''creates custom mapping settings, use mpl.rcdefaults() to default'''
    # can find out what is default by pprint(mpl.rcParams)
    mpl.rc('lines', linewidth=1.0, antialiased=True) # didn't change
    mpl.rc('patch', linewidth=0.5, facecolor='348ABD', edgecolor='eeeeee', \
        antialiased=True)
    mpl.rc('font', family='monospace', size=10.0) #, weight="bold")
    mpl.rc('axes', facecolor='f4f4f4', edgecolor='black', linewidth=2, \
        grid=False, titlesize='x-large', labelsize='x-large', \
        labelcolor='black', axisbelow=True)#, \
        #color_cycle='348ABD, 7A68A6, A60628, 467821, CF4457, 188487, E24A33')
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
        left=0.07, right=0.95, bottom=0.1, \
        top=0.95)


def create_cmap():
    #make a red-to-green colormap 
    colors = ['#980909', '#5C9809']
    cmap = LinearSegmentedColormap.from_list('my_map', colors, N = 50)
    return(cmap)


def create_plot(rc, lat_0, lon_0):
    #create a figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    #use custom matplotlib mappings defined above
    rc()
    
    #create a Basemap
    m = Basemap(resolution='l',projection='cyl',
    #m = Basemap(resolution='l',projection='cea',
        lat_ts=0,lat_0=lat_0,lon_0=lon_0, ax = ax)
    
    #create cmap
    cmap = create_cmap()

    # Add Grid Lines
    m.drawparallels(np.arange(-90, 90, 10.), labels=[1,0,0,0], fontsize=8)
    merids = m.drawmeridians(np.arange(-180, 190, 10.), labels=[0,0,0,1], fontsize=8)
    
    # Add Coastlines, States, and Country Boundaries
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    
    # Add Title
    plt.title('SIF (OCO-2)')
    
    #Rotate the x labels
    for merid in merids.keys():
        try:
            merids[merid][1][0].set_rotation(45)
        except:
            pass

    return(fig, ax, m, cmap)
   

def plot_data(ax, cmap, lats, lons, data):
    #plot the footprint parallelograms that fall within the requested range of longitudes (from plt_x)
    #PLOT-COLOR APPROACH #1
    plot_data = []
    #PLOT-COLOR APPROACH #2
    patches = []
    colors = []
    
    for i in range(lons.shape[0]):
        f_p = Polygon(zip(lons[i,:], lats[i,:]))
        f_x, f_y = f_p.exterior.xy
        if (np.int32(f_x)<plt_x[1]).all() and (np.int32(f_x)>plt_x[0]).all():
            f_xi, f_yi = m(f_x, f_y)
    
            #PLOT-COLOR APPROACH #1
            polygon = mpl_Polygon(np.array(list(zip(f_xi, f_yi))), closed=True)
            patches.append(polygon)
            colors.append(cmap(SIF[i]))
    
            #PLOT-COLOR APPROACH #2
            plot_data.append(f_xi)
            plot_data.append(f_yi)
            plot_data.append(mpl.colors.rgb2hex(cmap(SIF[i])))
            #ax.fill(f_xi, f_yi, mpl.colors.rgb2hex(cmap(SIF[i])))
    
    #PLOT-COLOR APPROACH #1
    #ax.fill(*plot_data);
    
    #PLOT-COLOR APPROACH #2
    collection = PatchCollection(patches)
    ax.add_collection(collection)
    collection.set_color(colors)


#call all of the above functions to create my plot
files, dates = get_netcdf_files(data_dir)
start, stop = list(map(convert_iso_date, [start, stop]))
print('\nFOUND THE FOLLOWING FILES:\n%s' % '\n'.join(['\t'+x for n, x in enumerate(files) if start <= dates[n] <= stop]))
print('----------------------------------------------------------\n\n\n')
for n, f in enumerate(files):
    if n == 0:
        fig, ax, m, cmap = create_plot(rc=rc, lat_0=lat_0, lon_0=lon_0)
    #only plot data for the correct year (for now)
    if start <= dates[n] <= stop:
        print('\nPlotting date %s' % dates[n])
        lons, lats, SIF, lon_0, lat_0 = read_netcdf(f, null_val = -999999.0)
        plot_data(ax=ax, cmap=cmap, lats=lats, lons=lons, data=SIF)

#show the plot
plt.show()


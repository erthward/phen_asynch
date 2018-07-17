#!/usr/bin/python
#prelim_sif_seasonality_analysis.py

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as mpl_Polygon
from matplotlib.colors import LinearSegmentedColormap

#set a range of longitudes within which to plot any data
plt_x = (-129, -34)

#read in the netCDF file as a Dataset object
netcdf = './data/oco2_LtSIF_180101_B8100r_180221185654s.nc4'
fh = Dataset(netcdf, mode = 'r')

#get the longitudes and  latitudes (which are 4-tuples of the vertices of each footprint parallelogram)
#and SIF values (for only 757nm in this case)
lons = fh.variables['footprint_vertex_longitude'][:]
lats = fh.variables['footprint_vertex_latitude'][:]
#SIF = fh.variables['SIF_757nm'][:]
SIF = fh.variables['SIF_771nm'][:]

#set the null value to that used in the netCDF
null_val = -999999.0

#get lon_0 and lat_0 values to center the plot on
lon_0 = mean(plt_x)
lat_0 = lats[np.where(lats != null_val)].mean()

#create a figure
fig = plt.figure()
ax = fig.add_subplot(111)

#create a Basemap
m = Basemap(width=15000000,height=18000000,
    resolution='l',projection='stere',
    lat_ts=40,lat_0=lat_0,lon_0=lon_0)

#make a red-to-green colormap 
colors = ['#980909', '#5C9809']
cmap = LinearSegmentedColormap.from_list('my_map', colors, N = 50)

# Plot Data
#cs = m.pcolor(xi,yi,np.squeeze(SIF))

# Add Grid Lines
m.drawparallels(np.arange(-90, 90, 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-179.75, 179.75, 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Title
plt.title('SIF (OCO-2)')

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

# Add Colorbar
#cbar = m.colorbar(cs, location='bottom', pad="10%")
#cbar.set_label(NDVI_units)

#show the plot
plt.show()

#close the netCDF file
fh.close()



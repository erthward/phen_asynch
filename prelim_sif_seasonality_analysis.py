#!/usr/bin/python
#prelim_sif_seasonality_analysis.py

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon
from matplotlib.colors import LinearSegmentedColormap

plt_x = (-129, -119)
#plt_x = (-129, -34)


netcdf = './oco2_LtSIF_180101_B8100r_180221185654s.nc4'
fh = Dataset(netcdf, mode = 'r')

lons = fh.variables['footprint_vertex_longitude'][:]
lats = fh.variables['footprint_vertex_latitude'][:]
SIF = fh.variables['SIF_757nm'][:]

lon_0 = mean(plt_x)
#lon_0 = lons.mean()
lat_0 = 0
#lat_0 = lats.mean()

fig = plt.figure()
ax = fig.add_subplot(111)

m = Basemap(width=45000000,height=22000000,
    resolution='l',projection='stere',
    lat_ts=40,lat_0=lat_0,lon_0=lon_0)

#make cmap
colors = ['#980909', '#5C9809']
cmap = LinearSegmentedColormap.from_list('my_map', colors, N = 50)

# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.
#lon, lat = np.meshgrid(lons, lats)
#xi, yi = m(lon, lat)


# Plot Data
#cs = m.pcolor(xi,yi,np.squeeze(SIF))

# Add Grid Lines
m.drawparallels(np.arange(-90, 90, 10.), labels=[1,0,0,0], fontsize=10)
m.drawmeridians(np.arange(-179.75, 179.75, 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Colorbar
#cbar = m.colorbar(cs, location='bottom', pad="10%")
#cbar.set_label(NDVI_units)

# Add Title
plt.title('SIF (OCO-2)')

plt.show()

#plot the footprint parallelograms that fall between the min and max lon vals in plt_x
for i in range(lons.shape[0]):
    f_p = Polygon(zip(lons[i,:], lons[i,:]))
    f_x, f_y = f_p.exterior.xy
    if (np.int32(f_x)<plt_x[1]).all() and (np.int32(f_x)>plt_x[0]).all():
        f_xi, f_yi = m(f_x, f_y)
        ax.fill(f_x, f_y, mpl.colors.rgb2hex(cmap(SIF[i])))


#close the netCDF file
fh.close()



from geopy.distance import great_circle, geodesic
from math import radians, cos, sin, asin, sqrt
import numpy as np


# haversine formula
def haversine(p1, p2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # get lons and lats
    lat1, lon1 = p1
    lat2, lon2 = p2

    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371008 # Radius of earth in meters
    return c * r


# points
p1 = (-42.27, 105.57)
p2 = (-41.94, 103.29)
p3 = (-41.94, 102.29)
p4 = (-4.94, 103.29)

# GEE dists
d0 = 0.14461657069437356
d1 = 191633.7501130543
d2 = 273050.7660794817
d3 = 4157039.7502771514
ds = [d0, d1, d2, d3]

#BEFORE THE FIX...
old_d0 = 0
old_d1 = 190157.12246519406
old_d2 = 270905.7564674099
old_d3 = 4151442.4463638524

# get dists from p1
n = 0
for p in [p1, p2, p3, p4]:
    print()
    print('%i.:' % n)
    print('\tGEODESIC:', geodesic(p1, p).meters)
    print()
    print('\tGREAT CIRCLE:', great_circle(p1, p).meters)
    print()
    print('\tHAVERSINE:', haversine(p1, p))
    print()
    print('\tGEE:', ds[n], '\t(off by:', np.abs(ds[n] - haversine(p1, p)), 'm)')
    print('---------------------------------')
    n += 1

# compare to distances calculated in GEE:


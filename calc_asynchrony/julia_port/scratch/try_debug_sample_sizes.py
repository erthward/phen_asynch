import geopy
import numpy as np
import matplotlib.pyplot as plt
from math import radians, cos, sin, asin, sqrt

def hav_dist(foc_x, foc_y, coord_pair):
    other_x = coord_pair[0]
    other_y = coord_pair[1]
    foc_x, foc_y, other_x, other_y = map(radians, [foc_x, foc_y, other_x,
                                                   other_y])
    dlon = other_x - foc_x
    dlat = other_y - foc_y
    a = sin(dlat/2)**2 + cos(foc_y) * cos(other_y) * sin(dlon/2)**2
    c = 2*asin(sqrt(a))
    r = 6371008
    return c*r

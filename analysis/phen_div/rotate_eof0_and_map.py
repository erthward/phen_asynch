import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

eofs = rxr.open_rasterio('../../../results/maps/global_4_EOFs_coswts.tif')

# rescale each layer 0-1
for i in range(eofs.shape[0]):
    eofs[i] = (eofs[i]-eofs[i].min())/(eofs[i].max()-eofs[i].min())

# get X and Y meshgrids
X, Y = np.meshgrid(eofs.x, eofs.y)

# get array indicating pixels for which to flip sign
flip_me = Y<0

# get first EOF, then deepcopy and flip sign on all pixels south of equator
eof0 = eofs[0].values
eof0_flip = deepcopy(eof0)
eof0_flip[flip_me] = -eof0_flip[flip_me]

# create weighting arrays to use for weighted sum of eof0 and eof0_flip,
# where the original values will be used above equator, the flipped values will
# be used below -23.4394 lat (i.e., tropic of capricorn), and the use of the two
# will be cosine-weighted between 0 and -23.4394
max_flip_lat = 23.4394
min_flip_lat = -23.4394
if max_flip_lat is None:
    max_flip_lat = np.max(eofs.y.values)
if min_flip_lat is None:
    min_flip_lat = np.min(eofs.y.values)
wts_eof0 = np.float32(np.invert(flip_me))
wts_eof0_flip = np.float32(flip_me)
for i, y in enumerate(eofs.y.values):
    if y<max_flip_lat and y>= min_flip_lat:
        frac = (y - min_flip_lat)/(max_flip_lat - min_flip_lat)
        wt_eof0 = np.cos(frac * (np.pi/2))
        wts_eof0[i,:] = 1-wt_eof0
        wts_eof0_flip[i,:] = wt_eof0

# calculate the weighted-sum eof0 layer, then plot it
eof0_wt_sum = (eof0 * wts_eof0) + (eof0_flip * wts_eof0_flip)
eofs[0] = eof0_wt_sum

eofs.plot.imshow(col='band', col_wrap=2)
plt.show()


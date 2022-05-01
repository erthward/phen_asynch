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
eof0_flip[flip_me] = -eof0_flip[flip_me]+1

# create weighting arrays to use for weighted sum of eof0 and eof0_flip,
# where the original values will be used above equator, the flipped values will
# be used below -23.4394 lat (i.e., tropic of capricorn), and the use of the two
# will be cosine-weighted between 0 and -23.4394
#max_flip_lat = 23.4394
#min_flip_lat = -23.4394
#if max_flip_lat is None:
#    max_flip_lat = np.max(eofs.y.values)
#if min_flip_lat is None:
#    min_flip_lat = np.min(eofs.y.values)
#wts_eof0 = np.float32(np.invert(flip_me))
#wts_eof0_flip = np.float32(flip_me)
#for i, y in enumerate(eofs.y.values):
#    if y<max_flip_lat and y>= min_flip_lat:
#        frac = (y - min_flip_lat)/(max_flip_lat - min_flip_lat)
#        wt_eof0 = np.cos(frac * (np.pi/2))
#        wts_eof0[i,:] = 1-wt_eof0
#        wts_eof0_flip[i,:] = wt_eof0

# calculate the weighted-sum eof0 layer, then plot it
#eof0_wt_sum = (eof0 * wts_eof0) + (eof0_flip * wts_eof0_flip)
eofs_trans_1 = deepcopy(eofs)
eofs_trans_1[0] = eof0_flip




# 04-30-2022: ANOTHER TAKE
eofs_trans_2 = deepcopy(eofs)
trans_arrays = []
# for 3 broad latitudinal zones (Americas, Europe and Africa, Oceania and Asia)
# calculate sliding-window latitudinal medians, then transform EOF 0 by
# calculating the absolute distance between each value and its respective
# latitudinal median
lon_slices = [[-180, -25], [-25, 61.5], [61.5, 180]]
for lon_slice in lon_slices:
    meds = eofs[0].sel(x=slice(*lon_slice)).median(dim='x').rolling(y=50).median()
    trans_arrays.append(np.abs(1-np.abs(eofs[0].sel(x=slice(*lon_slice)) - meds)))
eofs_trans_2[0] = np.hstack(trans_arrays)



# TAKE 3


# create a vertical vector that's 1 outside tropics, then goes from 1 to 0 and
# back, by cosine, within tropics
minmax_scale = lambda vals: (vals-np.min(vals))/(np.max(vals)-np.min(vals))
ys = eofs.y.values
wts = eofs.y*0
wts[ys<0] = 1
lats_in_tropics = eofs.sel(y=slice(23.4934, -23.4934)).y.values
# DO COSINE WEIGHTING
wts_in_tropics = 1-minmax_scale(np.cos(lats_in_tropics/(360)*(2*np.pi)))
# DO SIGMOID WEIGHTING
wts_in_tropics = 1/(1 + np.exp(-np.linspace(-10, 10, len(lats_in_tropics))))
wts[(ys >= -23.4934) * (ys < 23.4934)] = wts_in_tropics

# use that to weighted-sum rotated eof0 array and unrotated array
eofs_wt_sum = deepcopy(eofs)
eofs_wt_sum[0] = (wts*(1-eofs[0])) + ((1-wts)*eofs[0])



fig, axs = plt.subplots(1,2)
eofs_trans_1[:3].plot.imshow(ax=axs[0])
axs[0].set_title('RGB EOFs: weighted rotation')
eofs_wt_sum[:3].plot.imshow(ax=axs[1])
#eofs_wt_sum[[1,0,2]].plot.imshow(ax=axs[1])
#eofs_wt_sum[[1,2,0]].plot.imshow(ax=axs[1])
axs[1].set_title('RGB EOFs: tropical-sigmoid-weighted rotation')

fig.show()



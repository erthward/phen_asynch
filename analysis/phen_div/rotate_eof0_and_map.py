import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import geopandas as gpd
import palettable
import xycmap
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.colors import LinearSegmentedColormap
from colorsys import hls_to_rgb
from sklearn.metrics.pairwise import pairwise_distances_argmin
from sklearn.cluster import MiniBatchKMeans, KMeans
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from colormap import rgb2hex
from copy import deepcopy
import matplotlib as mpl
import rioxarray as rxr
from scipy import stats
import seaborn as sns
import pandas as pd
import numpy as np
#import xycmap
import os, re

# local imports
import helper_fns as hf


# TODO:
    # can I empirically figure out the amount of color 'deflation' as a
        # function of latitude, then multiply weighted-sum values by that to
        # 'reinflate'?
    # produce scree plot for EOF?
    # need to look into rotated EOFs as well?
    # need to center data (i.e., calc 'anomaly') prior to EOF, to make sure means = 0?

# load country boundaries
countries = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

# load the EOFs
eofs = rxr.open_rasterio(('/home/deth/Desktop/CAL/research/projects/'
                          'seasonality/results/maps/global_4_EOFs_coswts.tif'))

# rescale each layer 0-1
for i in range(eofs.shape[0]):
    eofs[i] = (eofs[i]-eofs[i].min())/(eofs[i].max()-eofs[i].min())

# create a vertical vector that's 0 in the far north, 1 in the far south,
# and progresses from 0 to 1 as a sigmoid function across the tropics
# back, by cosine, within tropics
minmax_scale = lambda vals: (vals-np.min(vals))/(np.max(vals)-np.min(vals))
ys = eofs.y.values
wts = eofs.y*0
wts[ys<0] = 1
lats_in_tropics = eofs.sel(y=slice(23.4934, -23.4934)).y.values
# create sigmoid weighting
# NOTE: MINMAX VAL CHOSEN HEURISTICALLY, TO MINIMIZE NOTICEABLE COLOR-WARPING
#       ARTEFACTS IN EQUATORIAL REGION 
minmax_val = 10
wts_in_tropics = minmax_scale(1/(1 + np.exp(-np.linspace(-minmax_val,
                                                         minmax_val,
                                                         len(lats_in_tropics)))))
wts[(ys >= -23.4934) * (ys < 23.4934)] = wts_in_tropics

# use that to weighted-sum the regular eof0 array and the inverted eof0 array
eofs_wt_sum = deepcopy(eofs)
eofs_wt_sum[0] = (wts*(1-eofs[0])) + ((1-wts)*eofs[0])

#fig, axs = plt.subplots(1,2)
#eofs[:3].plot.imshow(ax=axs[0])
#axs[0].set_title('RGB EOFs: no transform')
#eofs_wt_sum[:3].plot.imshow(ax=axs[1])
#eofs_wt_sum[[1,0,2]].plot.imshow(ax=axs[1])
#eofs_wt_sum[[1,2,0]].plot.imshow(ax=axs[1])
#axs[1].set_title('RGB EOFs: tropical sigmoid-weighted inversion')
#fig.show()

# write to raster file
#eofs_wt_sum = eofs_wt_sum.rio.write_crs(4326)
#eofs_wt_sum.rio.to_raster('./transformed_scaled_EOFs.tif')


# create bivariate colormap within luminance-difference on each of the two
# axes proportional to the percent variation explained by the two axes' EOFs
#pct_var_eof0 = 0.779
#pct_var_eof1 = 0.115
## define colors
#eof0_hue = 0.8
#eof0_saturation = 0.7
#eof1_hue = 0.3
#eof1_saturation = 0.7
#eof0_col0 = hls_to_rgb(eof0_hue, 0.2, eof0_saturation)
#eof0_col1 = hls_to_rgb(eof0_hue, 0.9, eof0_saturation)
#eof1_col0 = hls_to_rgb(eof1_hue, 0.2, eof1_saturation)
#eof1_col1 = hls_to_rgb(eof1_hue, 0.2+((pct_var_eof1/pct_var_eof0)*0.7), eof1_saturation)
#f, a = plt.subplots(1,1)
#a.scatter([0,1,0,1], [0,0,1,1], c=[eof0_col0, eof0_col1, eof1_col0, eof1_col1])
## NOTE: chosen manually
#ur_corner_col = '#a367c2'
#
#
## bivariate colormap
##fig2, ax = plt.subplots(1,1, figsize=(16,8))
##n = (50, 50)
##xcolors = ['#fae034', '#dcfa34', '#abfa34', '#52fa34', '#34fa9a', '#34fae0',
##           '#34cffa', '#3480fa', '#3445fa', '#6f34fa']
##xcmap = LinearSegmentedColormap.from_list('custom', xcolors)
##ycolors = ['#000000', '#ffffff']
##ycmap = LinearSegmentedColormap.from_list('custom', ycolors)
##cmap = xycmap.mean_xycmap(xcmap=xcmap, ycmap=ycmap, n=n)
#
## either pretty yellow, blue, and maroon colormap...
#corner_colors = ("#d1d1d1", "#f5e102", "#0098d9", "#87054d")
#
## or a data-driven colormap based on the code above...
##corner_colors = (eof0_col0, eof0_col1, eof1_col1, ur_corner_col)
#
## or a data-driven colormap hand-chosen with Google's hex color picker
## (with logic displayed in HSL values and arithmetic on luminosity value,
##  in comments to right)
##corner_colors = ("#014a46", # HSL = 177deg, 97%, 15%
##                 "#9df1fa", # HSL = 186deg, 90%, 80%
##                 "#7d7102", # HSL = 54deg, 96%, 25% (25% = 15% + (.115/.779)*65%)
##                 "#9ff9c4", # HSL = 145deg, 88%, 80%,
##                )
#cmap = xycmap.custom_xycmap(corner_colors=corner_colors, n=n)
#
#eof0_vals = eofs_wt_sum[0].values.ravel()
#eof1_vals = eofs_wt_sum[1].values.ravel()
#non_nans = np.where(np.invert(np.isnan(eof0_vals)) *
#                    np.invert(np.isnan(eof1_vals)))[0]
#sx=eof0_vals[np.invert(np.isnan(eof0_vals))]
#sy=eof1_vals[np.invert(np.isnan(eof1_vals))]
#colors = xycmap.bivariate_color(sx=sx, sy=sy, cmap=cmap)
#colors_arr = np.ones([len(eof0_vals), 3]) * np.nan
#colors_arr[non_nans] = np.array([c[:3] for c in colors])
#colors_arr = colors_arr.reshape([*eofs_wt_sum[0].shape, 3])
#colors_xr = deepcopy(eofs_wt_sum[:3])
#colors_xr[:,:,:] = np.swapaxes(np.swapaxes(colors_arr, 1, 2), 0, 1)
#
#countries.to_crs(4326).plot(color='none',
#                            linewidth=0.5,
#                            edgecolor='#4a4a4a',
#                            ax=ax,
#                            zorder=0,
#                           )
#colors_xr.plot.imshow(ax=ax)
#ax.set_xticks(())
#ax.set_yticks(())
#ax.set_xlabel('')
#ax.set_ylabel('')
#ax.set_title('')

# TODO: HOW TO ROTATE BIVAR LEGEND 45 DEG?
#plot_extents = -160, -110, -55, -20
#transform = Affine2D().rotate_deg(45)
#helper = floating_axes.GridHelperCurveLinear(transform, plot_extents)
#cax = floating_axes.FloatingSubplot(fig2, 111, grid_helper=helper)

#n_lgd = (8,8)
#cax = fig2.add_axes([0.07, 0.07, 0.12, 0.24])
#lgd_vals = np.meshgrid(*[np.linspace(0, 1, i) for i in n_lgd])
#lgd_cols = xycmap.bivariate_color(*[grid.ravel() for grid in lgd_vals], cmap=cmap)
#lgd_cols_arr = np.array([c[:3] for c in lgd_cols]).reshape([*n_lgd, 3])
#cax.imshow(lgd_cols_arr)
#cax.invert_yaxis()
#cax.set_xticks(np.linspace(0, n_lgd[0]-1, 5), np.round(np.linspace(0, 1, 5), 1))
#cax.set_yticks(np.linspace(0, n_lgd[1]-1, 5), np.round(np.linspace(0, 1, 5), 1))
#cax.set_xlabel('EOF 1', fontdict={'fontsize': 16})
#cax.set_ylabel('EOF 2', fontdict={'fontsize': 16})
#cax.tick_params(labelsize=9)
#fig2.subplots_adjust(left=0, bottom=0, right=1, top=1)
#fig2.show()
#
#fig2.savefig('phen_bivar_plot.png', dpi=700)





######################
# MAP K-MEANS CLUSTERS

# load EOFs
vals = eofs_wt_sum[:3].values.swapaxes(0,1).swapaxes(1,2)

# num clusters and batch size for mini-batch K-means
K_max = 12
batch_size = 40

# melt all 3 layers into an Nx3 array
X = vals.reshape([np.product(vals.shape[:2]), 3])
non_nans = np.where(np.sum(np.isnan(X), axis=1) == 0)[0]
X_sub = X[non_nans, :]

# simple K-means clustering
# (but using mini-batch,because it's faster)
# with within-cluster SoS recorded for each K
wcss = []
for K in range(1, K_max+1):
    mbk = MiniBatchKMeans(
        init="k-means++",
        n_clusters=K,
        batch_size=batch_size,
        n_init=10,
        max_no_improvement=10,
        verbose=0,
     )
    mbk.fit(X_sub)
    wcss.append(mbk.inertia_)


# set 'optimal' K
# (not really the elbow, but just a number that should help interpret
#  the overall spread of colors across the 3d space)
fig, ax = plt.subplots(1,1)
ax.plot(range(1, K_max+1), wcss)
opt_K = 11
ax.plot([opt_K, opt_K], ax.get_ylim())


# run K means for that value
mbk = MiniBatchKMeans(
        init="k-means++",
        n_clusters=opt_K,
        batch_size=batch_size,
        n_init=10,
        max_no_improvement=10,
        verbose=0,
     )
mbk.fit(X_sub)

# get colors for cluster centers
center_colors = [rgb2hex(*np.int32(mbk.cluster_centers_[i,:]*255)) for i in range(opt_K)]

# map the clusters
cluster_arr = np.ones([X.shape[0], 1])*np.nan
cluster_arr[non_nans,0] = mbk.labels_
clusters = cluster_arr.reshape(vals[:,:,0].shape)
clusters_rxr = deepcopy(eofs[0])
clusters_rxr[:,:] = clusters
fig, ax = plt.subplots(1,1, figsize=(16,8))
cmap = mpl.colors.ListedColormap(center_colors)
countries.to_crs(4326).plot(color='none',
                            linewidth=1,
                            edgecolor='black',
                            alpha=0.7,
                            ax=ax,
                            zorder=1,
                           )
clusters_rxr.plot.imshow(ax=ax, cmap=cmap, alpha=1, zorder=0)
fig.show()



# take 1/100th of the data in the X matrix,
# to be scattered as smaller, transparent points around the cluster centers
X_4scat = X_sub[range(0, X_sub.shape[0], 100), ]
X_4scat_colors = [rgb2hex(*np.int32(X_4scat[i,:]*255)) for i in range(X_4scat.shape[0])]


# plot clusters in 3D space, colored by RGB,
# and make interactive (to be able to twirl)
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X_4scat[:,0], X_4scat[:,1], X_4scat[:,2],
           c=X_4scat_colors, edgecolor='black', s=15, alpha=0.2)
ax.scatter(mbk.cluster_centers_[:,0],
           mbk.cluster_centers_[:,1],
           mbk.cluster_centers_[:,2],
           c=center_colors, edgecolor='black', s=350)
ax.set_xlabel('EOF 1', fontdict={"fontsize": 16})
ax.set_ylabel('EOF 2', fontdict={"fontsize": 16})
ax.set_zlabel('EOF 3', fontdict={"fontsize": 16})


##########################################################
# plot characteristic phenological cycles for each cluster

# load the coeffs
coeffs = rxr.open_rasterio(('/home/deth/Desktop/CAL/research/projects/'
                            'seasonality/results/maps/NIRv_global_coeffs.tif'))

# reshape the coeffs raster equivalent to how I reshaped eofs raster above
coeff_vals = coeffs.values.swapaxes(0,1).swapaxes(1,2)
coeff_X = coeff_vals.reshape([np.product(vals.shape[:2]), 5])
coeff_X_sub = coeff_X[non_nans, :]

# create a wts array and also reshape that
wts_vals = np.array([wts]*eofs[0].shape[1]).T
wts_X = wts_vals.reshape([np.product(wts_vals.shape[:2]), 1])
wts_X_sub = wts_X[non_nans, :]

# get harmonic regression design matrix
dm = hf.make_design_matrix()

# create figure and map RGB map
fig = plt.figure(figsize=(20,8))
grid_dim = int(np.ceil(opt_K/2))
gs = fig.add_gridspec(grid_dim, grid_dim)
ax_map = fig.add_subplot(gs[:, 1:(grid_dim-1)])
cmap = mpl.colors.ListedColormap(center_colors)
countries.to_crs(4326).plot(color='none',
                            linewidth=1,
                            edgecolor='black',
                            alpha=0.7,
                            ax=ax_map,
                            zorder=1,
                           )
eofs_wt_sum[:3].plot.imshow(ax=ax_map, cmap=cmap, alpha=1, zorder=0)

# for each cluster center, find N cells closest to it in RGB space
N = 1000
for i, c in enumerate(mbk.cluster_centers_):
    ax_clust = fig.add_subplot(gs[i%grid_dim, 0+((grid_dim-1)*(i>=grid_dim))])
    dists = []
    clust_member_coords = X_sub[mbk.labels_ == i,:]
    clust_member_coeffs = coeff_X_sub[mbk.labels_ == i,:]
    clust_member_wts = wts_X_sub[mbk.labels_ == i,:]
    for coords in clust_member_coords:
        dist = hf.calc_euc_dist(c, coords)
        dists.append(dist)
    closest = np.argsort(dists)[:100]
    # for each of those cells, get fitted, standardized seasonal TS,
    # then plot it
    tss = []
    for ind in closest:
        member_coeffs = clust_member_coeffs[ind]
        member_wt = clust_member_wts[ind]
        ts = hf.calc_time_series(member_coeffs, dm)
        # standardize ts
        ts = hf.standardize_array(ts)
        # rotate the ts the fraction of 180 days
        # indicated by the latitude's sigmoid weight
        ts_rot = np.concatenate((ts[int(182.5*(member_wt)):],
                                 ts[:int(182.5*(member_wt))]))
        assert(len(ts_rot) == len(ts))
        ax_clust.plot(ts_rot, color=center_colors[i], alpha=0.01)
        tss.append(ts_rot)
    # plot median and 5th and 95th percentiles
    tss_arr = np.array(tss)
    ax_clust.plot(np.median(tss_arr, axis=0),
                  color=center_colors[i],
                  linewidth=2.5)
    ax_clust.plot(np.percentile(tss_arr, 5, axis=0),
                  color=center_colors[i],
                  linewidth=2.5,
                  linestyle=':')
    ax_clust.plot(np.percentile(tss_arr, 95, axis=0),
                  color=center_colors[i],
                  linewidth=2.5,
                  linestyle=':')

import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from copy import deepcopy
import geopandas as gpd
import palettable
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


#############################################################################
# TODO:

    # can I empirically figure out the amount of color 'deflation' as a
        # function of latitude, then multiply weighted-sum values by that to
        # 'reinflate'?

    # produce scree plot for EOF?

    # need to look into rotated EOFs as well?

    # need to center data (i.e., calc 'anomaly') prior to EOF, to make sure means = 0?

#############################################################################



####################
# LOAD AND PREP DATA
####################

# load country boundaries
countries = gpd.read_file(('/home/deth/Desktop/CAL/research/projects/seasonality/'
                           'results/maps/NewWorldFile_2020.shp')).to_crs(8857)

# load the coeffs
coeffs = rxr.open_rasterio(('/home/deth/Desktop/CAL/research/projects/'
        'seasonality/results/maps/NIRv_global_coeffs.tif'))

# load the EOFs
eofs = rxr.open_rasterio(('/home/deth/Desktop/CAL/research/projects/'
    'seasonality/results/maps/global_4_EOFs_coswts.tif'))[:3]

# EOF percentages of variance explained
eofs_pcts = [70, 17, 8]

# define focal region bounding bboxes
reg_bboxes = {
          'CA': [-1.06e7, 4.6e6, -0.985e7, 2.87e6], # CA Mediterranean and Baja
          'QU': [1.32e7, -1.38e6, 1.373e7, -2.55e6], # N. Queensland peninsula
          'AM': [-6.3e6, -2e5, -5.8e6, -5.5e5], # upper Amazon river basin
          'AD': [-5.47e6, 3e5, -4.6e6, -6.35e5], # NE Arc of Deforestation
          'FL': [-7.525e6, 3.53e6, -7.25e6, 3.13e6], # South Florida
          #'CH': [-6.6e6, -3.6e6, -4.75e6, -5.75e6], # Chile Mediterranean
          #'SA': [1.3e6, -3.35e6, 2.7e6, -3.8e6], # S. Africa Mediterranean
          #'ME': [-1.5e6, 5.65e6, 7e6, -2.5e6], # actual Mediterranean
          #'AU': [1e7, -3.25e6, 1.36e7, -5.375e6], # Australia Mediterranean
         }
# define focal region gridspec indices
# NOTE: [[i_min, i_max], [j_min, j_max]]
reg_gsinds = {'QU': [[60, 90], [175, 195]],
              'AM': [[40, 60], [10, 30]],
              'AD': [[60, 90], [35, 55]],
              'FL': [[30, 60], [160, 180]],
              'CA': [[0, 30], [0, 30]],
             }
reg_gsinds_lines = {'QU': [[90, 100], [175, 195]],
                    'AM': [[60, 70], [10, 30]],
                    'AD': [[90, 100], [35, 55]],
                    'FL': [[60, 70], [160, 180]],
                    'CA': [[30, 40], [5, 25]],
                   }
# NOTE: K VALUES WERE DETERMINED BY MANUAL INSPECTION OF SCREE PLOTS
#       USING THE run_kmeans_clust FN WITH scree=True
reg_K_vals = {'QU': 3,
          'AM': 4,
          'AD': 3,
          'FL': 3,
          'CA': 4,
         }

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

# get equal-area-projected EOFs rasters, for mapping
eofs_wt_sum_for_map = eofs_wt_sum.rio.write_crs(4326)
eofs_wt_sum_for_map = eofs_wt_sum_for_map.rio.reproject(8857)
# and truncate artificially inflated values along edges
# (not sure what's broken that's causing them)
for lyr in range(eofs_wt_sum_for_map.shape[0]):
    eofs_wt_sum_for_map[lyr] = eofs_wt_sum_for_map[lyr].where(
        eofs_wt_sum_for_map[lyr] < 2*eofs_wt_sum[lyr].max(), np.nan)

# get harmonic regression design matrix
dm = hf.make_design_matrix()



###########
# PLOT EOFs
###########

# plotting helper fns
def strip_axes(ax):
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('')


# define xlims for global map
global_xlim = (0.90 * eofs_wt_sum_for_map.x.min(),
               0.95 * eofs_wt_sum_for_map.x.max())


# create figure
dims = (20, 10)
fig_main = plt.figure(figsize=(20,10))
gs = fig_main.add_gridspec(*[10*dim for dim in dims[::-1]]) # NOTE: REV ORDER OF FIGSIZE
# maps EOFS 1 (lat-weighted sum), 2, and 3
eofs_rows = [[0, 30], [0, 15], [15, 30]]
eofs_cols = [[30, 100], [100, 170], [100,170]]
for i, rows, cols in zip(range(3), eofs_rows, eofs_cols):
    ax_eof = fig_main.add_subplot(gs[rows[0]:rows[1], cols[0]:cols[1]])
    if i == 0:
        eofs_for_map = eofs.rio.write_crs(4326)
        eofs_for_map = eofs_for_map.rio.reproject(8857)
        eofs_for_map[0] = eofs_for_map[0].where(
                eofs_for_map[0] < 2*eofs[0].max(), np.nan)
        eofs_for_map[i].plot.imshow(ax=ax_eof,
                                    cmap='coolwarm',
                                    add_colorbar=False,
                                    alpha=1,
                                    zorder=0,
                                   )
        del eofs_for_map
    else:
        eofs_wt_sum_for_map[i].plot.imshow(ax=ax_eof,
                                           cmap='coolwarm',
                                           add_colorbar=False,
                                           alpha=1,
                                           zorder=0,
                                          )
    countries.plot(color='none',
                   linewidth=0.5,
                   edgecolor='black',
                   alpha=0.7,
                   ax=ax_eof,
                   zorder=1,
                  )
    strip_axes(ax_eof)
    ax_eof.set_xlim(global_xlim)
    ax_eof.set_ylim(eofs_wt_sum_for_map.rio.bounds()[1::2])
    ax_eof.text(0.92*ax_eof.get_xlim()[0], 0.92*ax_eof.get_ylim()[0],
                'EOF %i\n%0.1f%%' % (i, eofs_pcts[i]),
                fontdict={'fontsize': 14})
# plot the RGB map of EOFs 1-3 together
ax_rgb = fig_main.add_subplot(gs[35:, 35:165])
eofs_wt_sum_for_map.plot.imshow(ax=ax_rgb,
                                    add_colorbar=False,
                                    alpha=1,
                                    zorder=0,
                                   )
countries.plot(color='none',
               linewidth=0.5,
               edgecolor='black',
               alpha=0.7,
               ax=ax_rgb,
               zorder=1,
              )
strip_axes(ax_rgb)
ax_rgb.set_xlim(global_xlim)
ax_rgb.set_ylim(eofs_wt_sum_for_map.rio.bounds()[1::2])


####################
# PLOT FOCAL REGIONS
#################### 

def run_kmeans_clust(eofs_rast, coeffs_rast, reg,
                     k=None, k_max=12,
                     batch_size=40, n_init=10, max_no_improvement=10, seed=None,
                     n_clust_neighs=50,
                     plot_scree=False,
                     map_clusters=False, ax_lines=None):
    if seed is not None:
        np.random.seed(seed)
    if plot_scree:
        k = [*range(1, k_max+1)]
    else:
        assert k is not None, "if plot_scree==False then k must not be None"
        if not map_clusters:
            assert ax_lines is not None, ("axes must be provided if clustered "
                                    "time series line plot is to be generated")
    # prep EOFs vals for K-means clustering
    eofs_vals = eofs_rast.values.swapaxes(0,1).swapaxes(1,2)
    eofs_X = eofs_vals.reshape([np.product(eofs_vals.shape[:2]), 3])
    eofs_non_nans = np.where(np.sum(np.isnan(eofs_X), axis=1) == 0)[0]
    eofs_X_sub = eofs_X[eofs_non_nans, :]
    # loop over K vals and display scree plot, if requested
    if plot_scree:
        # simple K-means clustering
        # (but using mini-batch,because it's faster)
        # with within-cluster SoS recorded for each K
        wcss = []
        for k_val in k:
            mbk = MiniBatchKMeans(
                init="k-means++",
                n_clusters=k_val,
                batch_size=batch_size,
                n_init=n_init,
                max_no_improvement=max_no_improvement,
                verbose=0,
             )
            mbk.fit(eofs_X_sub)
            wcss.append(mbk.inertia_)
        fig_scree, ax_scree = plt.subplots(1,1)
        ax_scree.plot(range(1, k_max+1), wcss)
        ax_scree.set_title(reg, fontdict={'fontsize': 24})
        fig_scree.show()
        return
    # else just run clustering for the given value of k
    else:
        # run K means for that value
        mbk = MiniBatchKMeans(
                init="k-means++",
                n_clusters=k,
                batch_size=batch_size,
                n_init=n_init,
                max_no_improvement=max_no_improvement,
                verbose=0,
             )
        mbk.fit(eofs_X_sub)
    # get colors for cluster centers
    center_colors = [rgb2hex(*np.int32(mbk.cluster_centers_[i,
                                                   :]*255)) for i in range(k)]
    cmap = mpl.colors.ListedColormap(center_colors)
    # prep the coeffs vals the same way
    coeffs_vals = coeffs_rast.values.swapaxes(0,1).swapaxes(1,2)
    coeffs_X = coeffs_vals.reshape([np.product(coeffs_vals.shape[:2]), 5])
    coeffs_X_sub = coeffs_X[eofs_non_nans, :]
    # map the clusters, if requested
    if map_clusters:
        cluster_arr = np.ones([eofs_X.shape[0], 1])*np.nan
        cluster_arr[eofs_non_nans,0] = mbk.labels_
        clusters = cluster_arr.reshape(eofs_vals[:,:,0].shape)
        clusters_rxr = deepcopy(eofs_rast[0])
        clusters_rxr[:,:] = clusters
        fig_clust_map, ax_clust_map = plt.subplots(1,1)
        cmap = mpl.colors.ListedColormap(center_colors)
        countries.to_crs(4326).plot(color='none',
                                    linewidth=1,
                                    edgecolor='black',
                                    alpha=0.7,
                                    ax=ax_clust_map,
                                    zorder=1,
                                   )
        clusters_rxr.plot.imshow(ax=ax_clust_map,
                                 cmap=cmap, alpha=1, zorder=0)
        fig_clust_map.show()
        return
    # else, for each clust cent, find n_clust_neighs cells closest in RGB space
    for i, c in enumerate(mbk.cluster_centers_):
        # TODO: make separate axes for line plot
        dists = []
        clust_member_coords = eofs_X_sub[mbk.labels_ == i,:]
        clust_member_coeffs = coeffs_X_sub[mbk.labels_ == i,:]
        for coords in clust_member_coords:
            dist = hf.calc_euc_dist(c, coords)
            dists.append(dist)
        closest = np.argsort(dists)[:n_clust_neighs]
        # for each of those cells, get fitted, standardized seasonal TS,
        # then plot it
        tss = []
        for ind in closest:
            member_coeffs = clust_member_coeffs[ind]
            ts = hf.calc_time_series(member_coeffs, dm)
            # standardize ts
            ts = hf.standardize_array(ts)
            tss.append(ts)
        # plot median and 5th and 95th percentiles
        tss_arr = np.array(tss)
        ax_lines.plot(np.median(tss_arr, axis=0),
                      color=center_colors[i],
                      linewidth=2.5)
        ax_lines.plot(np.percentile(tss_arr, 5, axis=0),
                      color=center_colors[i],
                      linewidth=1,
                      linestyle=':')
        ax_lines.plot(np.percentile(tss_arr, 95, axis=0),
                      color=center_colors[i],
                      linewidth=1,
                      linestyle=':')


def plot_bbox_rectangle(bounds, ax):
    xmin = bounds[0]
    ymin = bounds[3]
    width = bounds[2] - xmin
    height = ymin - bounds[1]
    rect = Rectangle((xmin, ymin),
                     width, -height,
                     facecolor='none',
                     edgecolor='black',
                     linewidth=2)
    ax.add_patch(rect)


# try to help manage memory usage
del eofs
del eofs_wt_sum

for reg, bbox in reg_bboxes.items():
    print('now plotting region: %s' % reg)
    # plot the focal region
    gsi = reg_gsinds[reg]
    gsi_lines = reg_gsinds_lines[reg]
    #height_ratio = -(bbox[3]-bbox[1])/(bbox[2]-bbox[0])
    ax_reg = fig_main.add_subplot(gs[gsi[0][0]:gsi[0][1],
                                gsi[1][0]:gsi[1][1]])
    ax_lines = fig_main.add_subplot(gs[gsi_lines[0][0]:gsi_lines[0][1],
                                  gsi_lines[1][0]:gsi_lines[1][1]])
    eofs_wt_sum_for_map.sel(x=slice(bbox[0], bbox[2]),
                            y=slice(bbox[1], bbox[3])).plot.imshow(ax=ax_reg,
                                                                   zorder=0,
                                                            add_colorbar=False,
                                                                  )
    countries.plot(ax=ax_reg,
                   color='none',
                   edgecolor='black',
                   linewidth=0.5,
                   alpha=0.8,
                   zorder=1,
                  )
    strip_axes(ax_reg)
    ax_reg.set_xlim(bbox[0], bbox[2])
    ax_reg.set_ylim(bbox[3], bbox[1])

    # subset the focal region's data and run clustering
    eofs_wt_sum_foc = eofs_wt_sum_for_map.sel(x=slice(bbox[0], bbox[2]),
                                              y=slice(bbox[1], bbox[3]),
                                             ).rio.reproject(4326)
    eofs_wt_sum_foc = eofs_wt_sum_foc.where(eofs_wt_sum_foc !=
                                            eofs_wt_sum_foc._FillValue, np.nan)
    coeffs_foc = coeffs.rio.reproject(8857).sel(x=slice(bbox[0], bbox[2]),
                                                y=slice(bbox[1], bbox[3]),
                                             ).rio.reproject(4326)
    coeffs_foc = coeffs_foc.where(coeffs_foc != coeffs_foc._FillValue, np.nan)
    # run K-means clustering
    K = reg_K_vals[reg]
    run_kmeans_clust(eofs_wt_sum_foc, coeffs_foc, reg,
                     k=K,
                     ax_lines=ax_lines,
                     batch_size=40,
                     seed=123456)
    # adjust focal-region plotting and add bounding boxes
    strip_axes(ax_lines)
    ax_lines.set_xticks(np.linspace(0, 365, 5),
                        ['Jan', 'Apr', 'Jul', 'Oct', 'Jan'])
    for month in np.linspace(0, 365, 13):
        ax_lines.axvline(month, color='black',
                         linewidth=0.25, linestyle=':', alpha=0.75, zorder=0)
    ax_lines.tick_params(labelsize=14, rotation=0)
    for axis in ['top','bottom','left','right']:
        ax_lines.spines[axis].set_linewidth(2)
    plot_bbox_rectangle(bbox, ax_rgb)

# show and save figure
fig_main.show()
fig_main.subplots_adjust(left=0.02, right=0.98, bottom=0.04, top=0.98)
fig_main.savefig('ch2_fig2_EOF_maps.png', dpi=700)

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
from sklearn.cluster import MiniBatchKMeans, KMeans, DBSCAN
from sklearn.mixture import GaussianMixture
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from colormap import rgb2hex
from copy import deepcopy
import matplotlib as mpl
import rioxarray as rxr
import xarray as xr
from scipy import stats
from scipy.signal import argrelextrema
import seaborn as sns
import pandas as pd
import numpy as np
#import xycmap
import os, re

# local imports
import helper_fns as hf


#############################################################################
# TODO:

    # think through and decide best clustering algo for this (Spectral? DBSCAN?)

    # can I empirically figure out the amount of color 'deflation' as a
        # function of latitude, then multiply weighted-sum values by that to
        # 'reinflate'?

    # produce scree plot for EOF?

    # need to look into rotated EOFs as well?

    # need to center data (i.e., calc 'anomaly') prior to EOF, to make sure means = 0?

#############################################################################


########
# PARAMS
########

# save figs?
save_it = False

# which dataset to use?
dataset = 'NIRv'
#dataset = 'SIF'

# which mask mode to use?
masking_mode = 'default'
#masking_mode = 'strict'
mask_filename_ext_dict = {'strict': '_STRICTMASK',
                          'default': ''}
mask_filename_ext = mask_filename_ext_dict[masking_mode]


# which clustering algo to use?
clust_algo = 'kmeans'

# create the plots to aid EOF interpretation?
run_eof_interpretation = False



# helpful viz fn
def compare_wt_nowt_rasts(nowt, wt, bands=[0,1,2]):
    assert len(bands) in [1, 3]
    if len(bands) == 1:
        bands = bands[0]
    fig, axs = plt.subplots(1,2)
    try:
        nowt[bands,:,:].plot.imshow(ax=axs[0])
    except AttributeError:
        pass
    try:
        wt[bands,:,:].plot.imshow(ax=axs[1])
    except AttributeError:
        pass
    return fig



####################
# LOAD AND PREP DATA
####################

# load country boundaries
countries = gpd.read_file(('/home/deth/Desktop/CAL/research/projects/seasonality/'
                           'results/maps/NewWorldFile_2020.shp')).to_crs(8857)

# load ITCZ shapefile
# NOTE: digitized from Li and Zeng 2005, as reproduced in Zhisheng et al. 2015
itcz = gpd.read_file(('/home/deth/Desktop/CAL/research/projects/seasonality/'
                      'seasonal_asynchrony/analysis/'
                      'ITCZ_li_zeng_2005_digitized.shp'))

# load the coeffs
coeffs = rxr.open_rasterio(('/home/deth/Desktop/CAL/research/projects/'
        'seasonality/results/maps/%s_global_coeffs%s.tif') % (dataset,
                                                            mask_filename_ext))


# load the EOFs
eofs = rxr.open_rasterio(('/home/deth/Desktop/CAL/research/projects/'
    'seasonality/results/maps/%s_global_4_EOFs_coswts%s.tif') % (dataset,
                                                        mask_filename_ext))[:3]

# EOF percentages of variance explained
eofs_pcts = [70, 17, 8]

# define focal region bounding bboxes
reg_bboxes = {
              'Qu': [1.32e7, -1.38e6, 1.373e7, -2.55e6], # N. Queensland peninsula
              'AD': [-5.47e6, 3.1e5, -4.6e6, -0.6e6], # NE Arc of Deforestation
              'CA': [-1.06e7, 4.6e6, -0.985e7, 2.87e6], # CA Mediterranean and Baja
              'SAm':[-7.74e6, -0.6e6, -3.34e6, -2.3e6], # S. America cross-section
              'Af': [0.95e6, -0.5e6, 4.8e6, -2.3e6], # African cross-section
              'FL': [-7.525e6, 3.53e6, -7.25e6, 3.13e6], # South Florida
              'SAf':[1.52e6, -3.8e6, 2.14e6, -4.325e6], # S. Africa Mediterranean
              'Au': [1.008e7, -3.86e6, 1.09e7, -4.4e6], # Australia Mediterranean
             }
# define focal region gridspec indices
# NOTE: [[i_min, i_max], [j_min, j_max]]
reg_gsinds = {
              'Qu': [[0, 30], [180, 200]],
              'AD': [[25, 45], [30, 50]],
              'CA': [[0, 40], [0, 30]],
              'SAm':[[65, 85], [22, 94]],
              'Af': [[65, 85], [84, 148]],
              'FL': [[60, 90], [0, 20]],
              'SAf':[[70, 90], [150, 170]],
              'Au': [[60, 80], [175, 195]],
             }
# define gridspec indices for axes to plot line plots
# NOTE: should always be 10 indices tall and 20 wide
reg_gsinds_lines = {
              'Qu': [[30, 40], [180, 200]],
              'AD': [[45, 55], [30, 50]],
              'CA': [[40, 50], [5, 25]],
              'SAm':[[85, 95], [48, 68]],
              'Af': [[85, 95], [106, 126]],
              'FL': [[90, 100], [0, 20]],
              'SAf':[[90, 100], [150, 170]],
              'Au': [[80, 90], [175, 195]],
                   }
# NOTE: K VALUES WERE DETERMINED BY MANUAL INSPECTION OF SCREE PLOTS
#       USING THE run_clust_analysis FN WITH scree=True
reg_K_vals = {
              'Qu': 3,
              'AD': 3,
              'CA': 4,
              'SAm':5,
              'Af': 5,
              'FL': 3,
              'SAf':5,
              'Au': 5,
             }
reg_letters = {
              'Qu': 'E.',
              'AD': 'D.',
              'CA': 'B.',
              'SAm':'H.',
              'Af': 'I.',
              'FL': 'C.',
              'SAf':'G.',
              'Au': 'F.',
            }
# locations of region letter labels, expressed in fractions along x,y axes from top left
reg_lett_locs = {
              'Qu': (0.05, 0.9),
              'AD': (0.8, 0.8),
              'CA': (0.05, 0.05),
              'SAm':(0.05, 0.05),
              'Af': (0.9, 0.87),
              'FL': (0.05, 0.05),
              'SAf':(0.02, 0.05),
              'Au': (0.05, 0.9),
            }
assert (len(reg_bboxes) == len(reg_gsinds) ==
        len(reg_gsinds_lines) == len(reg_K_vals))
assert ([*reg_bboxes.keys()] == [*reg_gsinds.keys()] ==
        [*reg_gsinds_lines.keys()] == [*reg_K_vals.keys()])


def minmax_scale(vals, min_out=0, max_out=1):
    assert (min_out >=0) and (max_out <=1)
    scaled = (vals-np.min(vals))/(np.max(vals)-np.min(vals))
    scaled = scaled * (max_out - min_out)
    if min_out > 0:
        scaled = min_out + scaled
    return scaled


# rescale each layer 0-1
for i in range(eofs.shape[0]):
    eofs[i] = (eofs[i]-eofs[i].min())/(eofs[i].max()-eofs[i].min())


# create array implementing weights across ITCZ:
# 1. get DataArrays containing y vals as data, x vals as coordinates
#    (to facilitate lookup of y corresponding to nearest x to each x in eofs)
itcz_das = []
for n in range(3):
    xs = [itcz.geometry[n].coords[i][0] for i in
          range(len(itcz.geometry[n].coords))]
    ys = [itcz.geometry[n].coords[i][1] for i in
          range(len(itcz.geometry[n].coords))]
    da = xr.DataArray(data=ys, dims=['x'], coords=dict(x=(['x'], xs)))
    itcz_das.append(da)

# create empty wts DataArray, and empty wts numpy array to later go into it
wts = eofs[0]*np.nan
wts_arr = np.ones(wts.shape) * np.nan
#wts_inflation = eofs[0]*np.nan
#wts_inflation_arr = np.ones(wts.shape) * np.nan
# loop across x coordinates in wts
use_mean_itcz = True
if use_mean_itcz:
    itcz_delta_lat = 5
for j, x in enumerate(wts.x):
    # get N and S ys bounding the ITCZ at this x
    if use_mean_itcz:
        y_corresponding = float(itcz_das[2].sel(x=x, method='nearest'))
        y_N = y_corresponding + itcz_delta_lat
        y_S = y_corresponding - itcz_delta_lat
    else:
        y_N = float(itcz_das[0].sel(x=x, method='nearest'))
        y_S = float(itcz_das[1].sel(x=x, method='nearest'))
    # get number of N and S pixels outside the ITCZ at this x; set to 1 and 0
    n_N = wts.sel(x=x, y=slice(np.max(wts.y), y_N)).size
    n_S = wts.sel(x=x, y=slice(y_S, np.min(wts.y))).size
    wts_arr[:n_N, j] = 1
    wts_arr[-n_S:, j] = 0
    # create linear interpolation between 0 and 1 across ITCZ
    interp = np.linspace(1, 0, len(wts.y) - n_N - n_S + 2)
    assert len(interp) == 2 + np.sum(np.isnan(wts_arr[:, j]))
    wts_arr[n_N-1:-n_S+1, j] = interp
    # determine maximum possible weighted-sum value at this x;
    # add reciprocal of that to inflation array
    #wt_sum_maxes = np.max(np.stack([n*wts_arr[:,j] +
    #    ((1-n)*(1-wts_arr[:,j])) for n in np.linspace(0, 1, 100)]), axis=0)
    #inflation_factors = 1/wt_sum_maxes
    #wts_inflation_arr[:, j] = inflation_factors
wts.values = wts_arr
#wts_inflation.values = wts_inflation_arr

# make the weighted-sum EOFs map
eofs_wt_sum = deepcopy(eofs)
for n in range(2):
    eofs_wt_sum[n] = ((wts*(eofs[n])) + ((1-wts)*(1-eofs[n])))# * wts_inflation

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
global_xlim = (0.80 * eofs_wt_sum_for_map.x.min(),
               0.95 * eofs_wt_sum_for_map.x.max())


# create EOF fig
fig_eof = plt.figure(figsize=(20,30))

# maps EOFS 1, 2, and 3 (all raw)
eofs_for_map = eofs.rio.write_crs(4326)
eofs_for_map = eofs_for_map.rio.reproject(8857)
for i in range(3):
    ax_eof = fig_eof.add_subplot(3,1,i+1)
    eofs_for_map[i] = eofs_for_map[i].where(
            eofs_for_map[i] < 2*eofs[i].max(), np.nan)
    eofs_for_map[i].plot.imshow(ax=ax_eof,
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
    ax_eof.set_ylim(eofs_for_map.rio.bounds()[1::2])
    ax_eof.text(0.92*ax_eof.get_xlim()[0], 0.92*ax_eof.get_ylim()[0],
                'EOF %i\n%0.1f%%' % (i+1, eofs_pcts[i]),
                fontdict={'fontsize': 38})
del eofs_for_map
fig_eof.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98,
                        hspace=0.05)
if save_it:
    fig_eof.savefig('ch3_fig_%s_EOF_maps%s.png' % (dataset,
                                                   mask_filename_ext), dpi=700)


# create main figure
dims = (20, 10)
fig_main = plt.figure(figsize=(20,10))
gs = fig_main.add_gridspec(*[10*dim for dim in dims[::-1]]) # NOTE: REV ORDER OF FIGSIZE


# plot the RGB map of EOFs 1-3 together
ax_rgb = fig_main.add_subplot(gs[:65, 35:165])
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
ax_rgb.text(ax_rgb.get_xlim()[0]+0.01*np.diff(ax_rgb.get_xlim())[0],
            ax_rgb.get_ylim()[0]+0.95*np.diff(ax_rgb.get_ylim())[0],
            'A.',
            size=20,
            weight='bold',
           )


####################
# PLOT FOCAL REGIONS
#################### 

def calc_euc_dist(vec1, vec2):
    dist = np.sqrt(np.sum([(vec1[i] - vec2[i])**2 for i in range(len(vec1))]))
    return dist


def get_ssdists(vals, cent):
    dists = np.array([calc_euc_dist(val, cent) for val in vals])
    return dists**2


def calc_inertia(labels, uniq_labels, data, centers):
    """
    NOTE: aping inertia as explained in the KMeans docs
          (i.e., sum of squared dists of samples from their clust centers)
    NOTE: NOT ACTUALLY FROM THEIR 'CLOSEST' CLUST CENTERS!
    """
    ssdists = []
    for label in uniq_labels:
        ssdist = get_ssdists(data[labels==label, :], centers[label, :])
        ssdists.append(ssdist)
    inertia = np.sum(np.concatenate(ssdists))
    return inertia


def run_kmeans_clust(data, K, batch_size, n_init, max_no_improvement):
    clust = MiniBatchKMeans(
        init="k-means++",
        n_clusters=K,
        batch_size=batch_size,
        n_init=n_init,
        max_no_improvement=max_no_improvement,
        verbose=0,
    )
    clust.fit(data)
    inertia = clust.inertia_
    centers = clust.cluster_centers_
    labels = clust.labels_
    return labels, centers, inertia


def run_dbscan_clust(data, eps, min_samples):
    clust = DBSCAN(eps=eps, min_samples=min_samples)
    labels = clust.fit_predict(data)
    uniq_labels = np.unique(labels)
    centers = np.stack([np.mean(data[labels == lab, :],
                                axis=0) for lab in uniq_labels])
    inertia = calc_inertia(labels, uniq_labels, data, centers)
    return labels, centers, inertia


def run_gaussian_mix_clust(data, n_clusters):
    model = GaussianMixture(n_components=n_clusters)
    model.fit(data)
    labels = model.predict(data)
    uniq_labels = np.unique(labels)
    centers = np.stack([np.mean(data[labels == lab, :],
                                axis=0) for lab in uniq_labels])
    inertia = calc_inertia(labels, uniq_labels, data, centers)
    return labels, centers, inertia


def run_clust_analysis(eofs_rast, coeffs_rast, reg, clust_algo,
              k=None, k_max=12,
              batch_size=40, n_init=10, max_no_improvement=10, seed=None,
              eps=1, min_samples=10,
              n_clust_neighs=50,
              plot_envelope=False,
              plot_scree=False,
              map_clusters=False, ax_lines=None):
    assert clust_algo in ['kmeans', 'gaussian_mix', 'dbscan']
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
            if clust_algo == 'kmeans':
                labels, centers, inertia = run_kmeans_clust(eofs_X_sub,
                                                    k_val,
                                                    batch_size,
                                                    n_init,
                                                    max_no_improvement,
                                                   )
            elif clust_algo == 'gaussian_mix':
                labels, centers, inertia = run_gaussian_mix_clust(eofs_X_sub, k)
            elif clust_algo == 'dbscan':
                labels, centers, inertia = run_dbscan_clust(eofs_X_sub,
                                                            eps,
                                                            min_samples)
            wcss.append(inertia)
        fig_scree, ax_scree = plt.subplots(1,1)
        ax_scree.plot(range(1, k_max+1), wcss)
        ax_scree.set_title(reg, fontdict={'fontsize': 24})
        fig_scree.show()
        return
    # else just run clustering for the given value of k
    else:
        # run K means for that value
        if clust_algo == 'kmeans':
            labels, centers, inertia = run_kmeans_clust(eofs_X_sub,
                                                k,
                                                batch_size,
                                                n_init,
                                                max_no_improvement,
                                               )
        elif clust_algo == 'gaussian_mix':
            labels, centers, inertia = run_gaussian_mix_clust(eofs_X_sub, k)
        elif clust_algo == 'dbscan':
             labels, centers, inertia = run_dbscan_clust(eofs_X_sub,
                                                         eps,
                                                         min_samples)
    center_colors = [rgb2hex(*np.int32(centers[i, :]*255)) for i in range(k)]
    cmap = mpl.colors.ListedColormap(center_colors)
    # prep the coeffs vals the same way
    coeffs_vals = coeffs_rast.values.swapaxes(0,1).swapaxes(1,2)
    coeffs_X = coeffs_vals.reshape([np.product(coeffs_vals.shape[:2]), 5])
    coeffs_X_sub = coeffs_X[eofs_non_nans, :]
    # map the clusters, if requested
    if map_clusters:
        cluster_arr = np.ones([eofs_X.shape[0], 1])*np.nan
        cluster_arr[eofs_non_nans,0] = labels
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
    for i, c in enumerate(centers):
        # TODO: make separate axes for line plot
        dists = []
        clust_member_coords = eofs_X_sub[labels == i,:]
        clust_member_coeffs = coeffs_X_sub[labels == i,:]
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
        if plot_envelope:
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
#del eofs
#del eofs_wt_sum

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
    run_clust_analysis(eofs_wt_sum_foc, coeffs_foc, reg, clust_algo,
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
    # add letters to region maps
    ax_reg.text(ax_reg.get_xlim()[0]+reg_lett_locs[reg][0]*np.diff(ax_reg.get_xlim())[0],
                ax_reg.get_ylim()[0]+reg_lett_locs[reg][1]*np.diff(ax_reg.get_ylim())[0],
                reg_letters[reg],
                size=20,
                weight='bold',
               )

# show and save figure
fig_main.show()
fig_main.subplots_adjust(left=0.02, right=0.98, bottom=0.04, top=0.98)
if save_it:
    fig_main.savefig('ch3_fig_%s_RGB_EOF_map%s.png' % (dataset,
                                                    mask_filename_ext), dpi=700)



################
# INTERPRET EOFS
################
if run_eof_interpretation:

    def standardize(ts):
        return (ts-np.min(ts))/(np.max(ts)-np.min(ts))


    def calc_phase(ts):
        # TODO: ACTUALLY CALCULATE PHASE OF MAX FROM COEFFS
        return np.where(ts == np.max(ts))[0][0]


    def calc_concentration(ts):
        # TODO: BETTER METRIC HERE
        stand_ts = standardize(ts)
        entropy = -np.sum(ts * np.log(ts))
        return entropy


    def calc_modality(ts):
        stand_ts = standardize(ts)
        maxes = argrelextrema(stand_ts,
                              np.greater,
                              mode='wrap',
                              order=1,
                             )[0]
        if len(maxes) == 0:
            # NOTE: jitter all poitns tiny bit to avoid lack of extrema if fitted
            #       points symmetrically straddle the numerical maximum
            maxes = argrelextrema(stand_ts+np.random.normal(0, 0.000001, stand_ts.size),
                                  np.greater,
                                  mode='wrap',
                                  order=1,
                                 )[0]
        if len(maxes) not in [1, 2]:
            print(ts)
        assert len(maxes) in [1, 2], '%i maxes found!' % len(maxes)
        if len(maxes) == 1:
            return 0
        else:
            maxes_ratio = np.min(stand_ts[maxes])
            return maxes_ratio


    def calc_symmetry(ts):
        # TODO: FIGURE OUT A WAY TO CALCULATE THIS, THEN ADD TO DICT
        pass


    stats_fn_dict = {'phase': calc_phase,
                     'conc': calc_concentration,
                     'mod': calc_modality,
                     #'symm': calc_symmetry,
                    }

    eof_stats = {}
    #for i, eof in enumerate(eofs_wt_sum):
    for i, eof in enumerate(eofs):
        xs, ys = [arr.ravel() for arr in np.meshgrid(eof.x, eof.y)]
        vals = eof.values.ravel()
        lo_pctile = np.nanpercentile(vals, 1)
        hi_pctile = np.nanpercentile(vals, 99)
        lo_inds = np.where(vals<=lo_pctile)
        hi_inds = np.where(vals>=hi_pctile)
        lo_xs = xs[lo_inds]
        lo_ys = ys[lo_inds]
        hi_xs = xs[hi_inds]
        hi_ys = ys[hi_inds]
        lo_stats = {stat: [] for stat in stats_fn_dict.keys()}
        hi_stats = {stat: [] for stat in stats_fn_dict.keys()}
        for lo_x, lo_y in zip(lo_xs, lo_ys):
            lo_coeffs = coeffs.sel(x=lo_x, y=lo_y, method='nearest').values
            ts = hf.calc_time_series(lo_coeffs, dm)
            for stat, fn in stats_fn_dict.items():
                lo_stats[stat].append(fn(ts))
        for hi_x, hi_y in zip(hi_xs, hi_ys):
            hi_coeffs = coeffs.sel(x=hi_x, y=hi_y, method='nearest').values
            ts = hf.calc_time_series(hi_coeffs, dm)
            for stat, fn in stats_fn_dict.items():
                hi_stats[stat].append(fn(ts))
        lo_stats_df = pd.DataFrame(lo_stats)
        lo_stats_df['pctile'] = 'lo'
        hi_stats_df = pd.DataFrame(hi_stats)
        hi_stats_df['pctile'] = 'hi'
        stats_df = pd.concat((lo_stats_df, hi_stats_df))
        stats_df['eof'] = i + 1
        eof_stats[i] = stats_df
    df = pd.concat([*eof_stats.values()])


    fig_eof_interp = plt.figure(figsize=(10,10))
    gs = fig_eof_interp.add_gridspec(len(stats_fn_dict),len(eofs))
    for i, stat in enumerate(stats_fn_dict.keys()):
        for j in df['eof'].unique():
            ax = fig_eof_interp.add_subplot(gs[i, j-1])
            subdf = df[df['eof'] == j].loc[:, [stat, 'pctile']]
            sns.violinplot(x='pctile',
                           y=stat,
                           data=subdf,
                           ax=ax,
                          )
            stat_range = np.max(subdf[stat]) - np.min(subdf[stat])
            ax.set_ylim(np.min(subdf[stat]) - 0.1*stat_range,
                        np.max(subdf[stat]) + 0.1*stat_range)
            if (j-1) == 0:
                ax.set_ylabel(stat, fontdict={'fontsize': 22})
                ax.tick_params(labelsize=14)
            else:
                ax.set_ylabel('')
                ax.set_yticks(())
            if i == 0:
                ax.set_title('EOF %i' % j, fontdict={'fontsize': 22})
            else:
                ax.set_title('')
            if i == (len(stats_fn_dict)-1):
                ax.set_xlabel('pctile', fontdict={'fontsize': 16})
                ax.set_xticks([0, 1], ['lo', 'hi'])
                ax.tick_params(labelsize=14)
            else:
                ax.set_xlabel('')
                ax.set_xticks(())
    fig_eof_interp.subplots_adjust(bottom=0.08,
                                   top=0.92,
                                   left=0.1,
                                   right=0.98,
                                   hspace=0.05,
                                   wspace=0.05,
                                  )
    fig_eof_interp.show()
    if save_it:
        fig_eof_interp.savefig('ch3_%s_EOF_interpretation_fig%s.png' % (dataset,
                                                mask_filename_ext), dpi=600)

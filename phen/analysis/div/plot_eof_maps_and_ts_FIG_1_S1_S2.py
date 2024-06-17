import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patheffects as patheffects
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
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from colormap import rgb2hex
from copy import deepcopy
import matplotlib as mpl
import rioxarray as rxr
import xarray as xr
from scipy import stats
from scipy.signal import argrelmax
import seaborn as sns
import pandas as pd
import numpy as np
#import xycmap
import os, sys, re

# local imports
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf


########
# PARAMS
########

# which figure to plot?
#what_to_plot = 'fig_1'
what_to_plot = 'reg_figs'
#what_to_plot = 'fig_s3'
#what_to_plot = 'fig_s4'

# plotting params
partlabel_fontsize = 24

# write the weighted-sum EOF map to a raster file?
write_wt_sum_eofs_to_file = False

# save figs?
save_it = True

# which dataset to use?
dataset = 'NIRv'
#dataset = 'SIF'

# which mask mode to use?
masking_mode = 'default'
#masking_mode = 'strict'
mask_filename_ext_dict = {'strict': '_STRICTMASK',
                          'default': ''}
mask_filename_ext = mask_filename_ext_dict[masking_mode]

# use the 'folding' method to compare N and S hemispheres?
fold_it = True

# were ts normalized before EOF calculation?
normts = True
if normts:
    normts_file_substr = '_normts'
else:
    normts_file_substr = ''


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
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                       'NewWorldFile_2020.shp')).to_crs(8857)

# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR, f)).to_crs(8857))
subnational = pd.concat(subnational)

# load ITCZ shapefile
# NOTE: digitized from Li and Zeng 2005, as reproduced in Zhisheng et al. 2015
itcz = gpd.read_file(os.path.join(phf.DATA_DIR, 'itcz',
                                  'ITCZ_li_zeng_2005_digitized.shp'))

# load the coeffs
coeffs = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                        '%s%s_coeffs.tif') % (dataset,
                                                            mask_filename_ext))


# load the EOFs
eofs = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                      '%s_global_4_EOFs_sqrt_coswts%s%s.tif') %
                         (dataset, normts_file_substr, mask_filename_ext))[:3]

# EOF percentages of variance explained
# TODO: LOAD THESE AND TIME SERIES FROM CSV!
eofs_pcts = [np.nan]

# define focal region bounding bboxes
reg_bboxes = {
              'QLD': [1.32e7, -1.375e6, 1.373e7, -2.45e6], # far N. Queensland
              'Am': [-5.47e6, 3.1e5, -4.6e6, -0.6e6],     # mouth of Amazon
              'Ba': [-1.06e7, 4.5e6, -0.985e7, 2.87e6],   # Baja
              'GB': [-1.01e7, 5.5e6, -0.94e7, 4.6e6],     # Great Basin
              'Mad':[3.94e6, -1.4e6, 4.84e6, -3.3e6],     # Madagascar
              'Fl': [-7.525e6, 3.53e6, -7.25e6, 3.13e6],  # South Florida
              'SAf':[1.57e6, -4.0e6, 2.2e6, -4.34e6],     # S. Africa cape
              'Au': [1.009e7, -3.959e6, 1.092e7, -4.364e6],  # SW Australia Medit.
              'IT': [0.55e6, 5.62e6, 1.2e6, 5.35e6],      # Po Valley, Italy
              'CB': [-8.9e6, 5.95e6, -6.95e6, 3.6e6],     # US Corn Belt
             }


# NOTE: K VALUES WERE DETERMINED BY MANUAL INSPECTION OF SCREE PLOTS
#       USING THE run_clust_analysis FN WITH plot_scree=True
reg_K_vals = {
              'QLD': 3,
              'Am': 3,
              'Ba': 4,
              'GB': 3,
              'Mad':3,
              'Fl': 3,
              'SAf':4,
              'Au': 4,
              'IT': 3,
              'CB': 3,
             }


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

if fold_it:
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
    # NOTE: folding all 3 EOFs, on the assumption that a hemispheric signal is
    #       embedded in each, because the signal in EOFs 2 and 3 is only
    #       obvious when the RGB composite is displayed, but not in their
    #       individual maps (unlike EOF 1)
    for n in range(3):
        eofs_wt_sum[n] = ((wts*(eofs[n])) + ((1-wts)*(1-eofs[n])))# * wts_inflation

    # get equal-area-projected EOFs rasters, for mapping
    eofs_wt_sum_for_map = eofs_wt_sum.rio.write_crs(4326)
else:
    eofs_wt_sum_for_map = eofs_wt_sum.rio.write_crs(4326)
eofs_wt_sum_for_map = eofs_wt_sum_for_map.rio.reproject(8857)
# and truncate artificially inflated values along edges
# (not sure what's broken that's causing them)
for lyr in range(eofs_wt_sum_for_map.shape[0]):
    eofs_wt_sum_for_map[lyr] = eofs_wt_sum_for_map[lyr].where(
        eofs_wt_sum_for_map[lyr] < 2*eofs_wt_sum[lyr].max(), np.nan)

# if requested, write mapping-prepped EOFS to file
if write_wt_sum_eofs_to_file:
    x_vals = eofs_wt_sum_for_map['x'].values
    y_vals = eofs_wt_sum_for_map['y'].values
    bands = eofs_wt_sum_for_map['band'].values-1
    out_bands = []
    for band in bands:
        out_bands.append(eofs_wt_sum_for_map[band].data)
    out_da = xr.DataArray(out_bands,
                          coords={'y':y_vals,
                                  'x':x_vals,
                                  'band': bands,
                                 },
                          dims=['band', 'y', 'x'],
                         )
    out_da.rio.to_raster(('/media/deth/SLAB/diss/3-phn/GEE_outputs/final/'
                          'NIRv_global_4_EOF_s_sqrt_coswts_normts_FORMAP.tif'))


# get harmonic regression design matrix
dm = phf.make_design_matrix()


############
## PLOT EOFs
############
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


if what_to_plot == 'fig_s3':
    # create EOF fig
    fig_eof = plt.figure(figsize=(20,30))

    # make maps of EOFS 1, 2, and 3 (all raw)
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
        subnational.plot(color='none',
                         linewidth=0.3,
                         edgecolor='black',
                         alpha=0.5,
                         ax=ax_eof,
                         zorder=1,
                        )
        countries.plot(color='none',
                       linewidth=0.5,
                       edgecolor='black',
                       alpha=0.7,
                       ax=ax_eof,
                       zorder=2,
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
        fig_eof.savefig('FIG_S1_%s_EOF_maps%s%s.png' % (dataset,
                    mask_filename_ext, ('_RAW'*(not fold_it))), dpi=600)



#############################
# PLOT UNTRANSFORMED RGB MAPS
#############################
if what_to_plot == 'fig_s4':
    # create untransformed RGB figure
    eofs_for_map = eofs.rio.write_crs(4326)
    eofs_for_map = eofs_for_map.rio.reproject(8857)
    fig_untrans = plt.figure(figsize=(20,20))
    ax_top = fig_untrans.add_subplot(2,1,1)
    ax_bot = fig_untrans.add_subplot(2,1,2)
    axs = {0: ax_top, 1: ax_bot}
    for row in range(2):
        ax = axs[row]
        subnational.plot(color='none',
                         linewidth=0.3,
                         edgecolor='black',
                         alpha=0.5,
                         ax=ax,
                         zorder=1,
                        )
        countries.plot(color='none',
                       linewidth=0.5,
                       edgecolor='black',
                       alpha=0.7,
                       ax=ax,
                       zorder=2,
                      )
        if row == 0:
            eofs_for_map.plot.imshow(ax=ax, zorder=0)
        elif row == 1:
            eofs_trans = deepcopy(eofs_for_map)
            # NOTE: folding all 3 EOFs, on the assumption that a hemispheric signal is
            #       embedded in each, rather than eyeballing whethere or not one is...
            for i in range(3):
                eofs_trans[i,:,:] = 1 - eofs_trans[i,:,:]
            eofs_trans = eofs_trans.where(np.abs(eofs_trans)<1e37)
            eofs_trans.plot.imshow(ax=ax, zorder=0)
        strip_axes(ax)
        ax.set_xlim(global_xlim)
        ax.set_ylim(eofs_for_map.rio.bounds()[1::2])
    fig_untrans.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98,
                                hspace=0.05)
    if save_it:
        fig_untrans.savefig('FIG_S2_untransformed_EOF_maps.png', dpi=600)

    del eofs_for_map



##############
# PLOT RGB MAP
##############


if what_to_plot == 'fig_1':

    # create main figure
    dims = (20, 10)
    fig_1 = plt.figure(figsize=(20,10))
    ax = fig_1.add_subplot(111)

    # plot the RGB map of EOFs 1-3 together
    eofs_wt_sum_for_map.plot.imshow(ax=ax,
                                        add_colorbar=False,
                                        alpha=1,
                                        zorder=0,
                                       )
    subnational.plot(color='none',
                     linewidth=0.3,
                     edgecolor='black',
                     alpha=0.5,
                     ax=ax,
                     zorder=1,
                    )
    countries.plot(color='none',
                   linewidth=0.5,
                   edgecolor='black',
                   alpha=0.7,
                   ax=ax,
                   zorder=2,
                  )

    strip_axes(ax)
    ax.set_xlim(global_xlim)
    ax.set_ylim(eofs_wt_sum_for_map.rio.bounds()[1::2])
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')

    # save figure
    fig_1.subplots_adjust(left=0.01,
                          right=0.99,
                          bottom=0.04,
                          top=0.98)
    if save_it:
        fig_1.savefig('FIG_1_%s_RGB_EOF_map%s%s.png' % (dataset,
                            mask_filename_ext, ('_RAW'*(not fold_it))), dpi=700)


#################
# PLOT FOCAL MAPS
#################


if what_to_plot == 'reg_figs':

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
        # NOTE: or, if this is the Great Basin, create the cluster raster,
        #       then overlay with the cheatgrass data and run a simple ANOVA
        #       and Tukey's HSD to test differences in modeled percent annual
        #       herbaceous cover between clusters
        if map_clusters or reg == 'GB':
            cluster_arr = np.ones([eofs_X.shape[0], 1])*np.nan
            cluster_arr[eofs_non_nans,0] = labels
            clusters = cluster_arr.reshape(eofs_vals[:,:,0].shape)
            clusters_rxr = deepcopy(eofs_rast[0])
            clusters_rxr[:,:] = clusters
            if map_clusters:
                fig_clust_map, ax_clust_map = plt.subplots(1,1)
                cmap = mpl.colors.ListedColormap(center_colors)
                subnational.plot(color='none',
                     linewidth=0.3,
                     edgecolor='black',
                     alpha=0.5,
                     ax=ax_clust_map,
                     zorder=1,
                    )
                countries.to_crs(4326).plot(color='none',
                                            linewidth=0.5,
                                            edgecolor='black',
                                            alpha=0.7,
                                            ax=ax_clust_map,
                                            zorder=2,
                                           )
                try:
                    clusters_rxr.plot.imshow(ax=ax_clust_map,
                                             cmap=cmap, alpha=1, zorder=0)
                # ignore that annoying rioxarray error
                # ("'tuple' object has no attribute 'startswith'")
                except AttributeError as e:
                    pass
                fig_clust_map.show()
            # load data on percent annual cover and run stats test of significant
            # difference between clusters
            if reg == 'GB':
                # data from: https://www.sciencebase.gov/catalog/item/5ec5159482ce476925eac3b7
                annuals = rxr.open_rasterio(('./WGA_Weighted_Mean_Annual_'
                                    'Herbaceous_Cover_2016_2017_2018_AGG5KM.tif'),
                                            masked=True)
                annuals = annuals.rio.reproject_match(clusters_rxr)
                cluster_vals = {i: annuals[0].values[clusters_rxr == i].ravel(
                                                ) for i in np.unique(clusters_rxr)}
                cluster_vals = {i: vals[pd.notnull(vals)] for i,
                                vals in cluster_vals.items() if pd.notnull(i)}
                # run ANOVA
                fval, pval = stats.f_oneway(*cluster_vals.values())
                # run Tukey HSD
                tuk = pairwise_tukeyhsd(
                    endog=np.array([val for val_list in cluster_vals.values(
                                                        ) for val in val_list]),
                    groups=np.array([clust for clust, val_list in cluster_vals.items(
                                                        ) for val in val_list]))
                print('\n\n==========\nGreat Basin cheatgrass statistical results:')
                print('\n\tANOVA: F-val: %0.2f; p-val: %0.2f\n\n' % (fval, pval))
                print('\tgroup means:\n\t%s' % '\t'.join(['%i: %0.2f' % (i,
                        np.mean(vals)) for i, vals in cluster_vals.items()]))
                print('\n\n\tTukey HSD: %s' % str(tuk))
                print('\n\ncluster colors:\n')
                for clust, color in zip(sorted(np.unique(clusters_rxr)),
                                        center_colors):
                    print(f'cluster {clust}: {color}\n')
                print('\n\n==========\n\n')
            if map_clusters:
                return
        # else, for each clust cent, find n_clust_neighs cells closest in RGB space
        for i, c in enumerate(centers):
            # TODO: make separate axes for line plot
            dists = []
            clust_member_coords = eofs_X_sub[labels == i,:]
            clust_member_coeffs = coeffs_X_sub[labels == i,:]
            for coords in clust_member_coords:
                dist = phf.calc_euc_dist(c, coords)
                dists.append(dist)
            closest = np.argsort(dists)[:n_clust_neighs]
            # for each of those cells, get fitted, standardized seasonal TS,
            # then plot it
            tss = []
            for ind in closest:
                member_coeffs = clust_member_coeffs[ind]
                ts = phf.calc_time_series(member_coeffs, dm)
                # standardize ts
                ts = phf.standardize_array(ts)
                tss.append(ts)
            # plot median and 5th and 95th percentiles
            tss_arr = np.array(tss)
            ax_lines.plot(np.median(tss_arr, axis=0),
                          color=center_colors[i],
                          linewidth=3.5)
            if plot_envelope:
                ax_lines.fill_between(x=np.linspace(*ax_lines.get_xlim(),
                                                  tss_arr.shape[1]),
                                      y1=np.percentile(tss_arr, 5, axis=0),
                                      y2=np.percentile(tss_arr, 95, axis=0),
                                      color=center_colors[i],
                                      alpha=0.7,
                                     )
                #ax_lines.plot(np.percentile(tss_arr, 5, axis=0),
                #              color=center_colors[i],
                #              linewidth=1,
                #              linestyle=':')
                #ax_lines.plot(np.percentile(tss_arr, 95, axis=0),
                #              color=center_colors[i],
                #              linewidth=1,
                #              linestyle=':')
        ax_lines.set_xlim((0, 365))
        return center_colors


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


    def add_phen_labs(ax, reg, text_size=12, mark_size=65, hspace_frac=0.22):
        """
        Add colored numeric labels to a phenology line plot
        """
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        n_lines = len(ax.get_lines())
        if n_lines <= 3:
            hspace_top_offset=1.05
        elif n_lines == 4:
            hspace_top_offset = 0.65
        relmaxes = [argrelmax(l.get_ydata(),
                              order=40,
                              mode='clip',
                             )[0] for l in ax.get_lines()]
        in_90pctle = [ax.get_lines()[i].get_ydata()[relmaxes[i]] ==
                      np.max(ax.get_lines()[i].get_ydata()) for
                      i in range(len(relmaxes))]
        relmaxes = [relmaxes[i][in_90pctle[i]][0] for i in range(len(relmaxes))]
        line_order = np.argsort(relmaxes)
        # NOTE: reordering line colors for South Africa so that color numbers
        #       are more intuitive in the figure
        if reg == 'SAf':
            line_order = [line_order[-1]] + list(line_order[:-1])
        colors = [ax.get_lines()[n].get_color() for n in line_order]
        for n, color in enumerate(colors):
            lab = n+1
            ax.text(400,
                    (ylim[1]-
                     (n*hspace_frac*np.diff(ylim))-
                     (0.5*hspace_frac)-
                     hspace_top_offset+0.1),
                    str(lab),
                    color='black',
                    size=text_size,
                    clip_on=False,
                   )
            ax.scatter(373,
                       (ylim[1]-
                        (n*hspace_frac*np.diff(ylim))-
                        hspace_top_offset),
                       c=color,
                       edgecolor='black',
                       marker='s',
                       s=mark_size,
                       clip_on=False)
        ax.set_xlim(xlim)


    def get_fig_dims_from_reg_bbox(bbox, fig_max_dim=7):
        bbox_xdim = bbox[2] - bbox[0]
        bbox_ydim = bbox[1] - bbox[3]
        bbox_maxdim = np.max([bbox_xdim, bbox_ydim])
        fig_xfact = bbox_xdim/bbox_maxdim
        fig_yfact = bbox_ydim/bbox_maxdim
        fig_dims = [fig_max_dim * fact for fact in [fig_xfact, fig_yfact]]
        return fig_dims

    for reg, bbox in reg_bboxes.items():

        print('now plotting region: %s' % reg)

        # create regional figure
        dims = get_fig_dims_from_reg_bbox(bbox)
        fig_reg, ax_reg = plt.subplots(1, figsize=dims)
        fig_lines, ax_lines = plt.subplots(1, figsize=(4,1))

        # plot the focal region
        eofs_wt_sum_for_map.sel(x=slice(bbox[0], bbox[2]),
                                y=slice(bbox[1], bbox[3])).plot.imshow(ax=ax_reg,
                                                                       zorder=0,
                                                                add_colorbar=False,
                                                                      )
        subnational.plot(color='none',
                         linewidth=0.3,
                         edgecolor='black',
                         alpha=0.5,
                         ax=ax_reg,
                         zorder=1,
                        )
        countries.plot(ax=ax_reg,
                       color='none',
                       edgecolor='black',
                       linewidth=0.5,
                       alpha=0.8,
                       zorder=2,
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
        clust_colors = run_clust_analysis(eofs_wt_sum_foc,
                                          coeffs_foc, reg,
                                          clust_algo,
                                          plot_envelope=False,
                                          k=K,
                                          ax_lines=ax_lines,
                                          batch_size=40,
                                          seed=123456)
        # add lineplot labels
        add_phen_labs(ax_lines, reg)
        # adjust focal-region plotting and add bounding boxes
        strip_axes(ax_lines)
        ax_lines.set_xticks(np.linspace(0, 365, 5))
        ax_lines.set_xticklabels(['Jan', 'Apr', 'Jul', 'Oct', 'Jan'])
        for month in np.linspace(0, 365, 13):
            ax_lines.axvline(month, color='black',
                             linewidth=0.25, linestyle=':', alpha=0.75, zorder=0)
        ax_lines.tick_params(labelsize=10, rotation=0)
        for axis in ['top','bottom','left','right']:
            ax_lines.spines[axis].set_linewidth(2)

        # save figures
        fig_reg.subplots_adjust(left=0.01,
                                 right=0.99,
                                 bottom=0.04,
                                 top=0.98)
        fig_lines.subplots_adjust(left=0.04,
                                 right=0.94,
                                 bottom=0.24,
                                 top=0.97)

        if save_it:
            fig_reg.savefig('FIG_REG_%s_%s_RGB_EOF_reg_map%s%s.png' % (reg, dataset,
                                mask_filename_ext, ('_RAW'*(not fold_it))), dpi=700)
            fig_lines.savefig('FIG_REG_%s_LINES_%s_RGB_EOF_reg_map%s%s.png' % (reg, dataset,
                                mask_filename_ext, ('_RAW'*(not fold_it))), dpi=700)

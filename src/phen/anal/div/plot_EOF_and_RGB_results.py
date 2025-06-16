import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patheffects as patheffects
from matplotlib.gridspec import GridSpec
from copy import deepcopy
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
import palettable
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from sklearn.metrics.pairwise import pairwise_distances_argmin
from sklearn.cluster import MiniBatchKMeans
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
import numpy as np
import os, sys, re
import pickle

# local imports
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf



########
# PARAMS
########

# what figures to plot?
try:
    what_to_plot = sys.argv[1]
    assert what_to_plot in ['main_rgb_map', 'reg_figs',
                            'eof_summ_fig', 'raw_rgb_maps']
except Exception:
    what_to_plot = input("what_to_plot: ")
assert what_to_plot in ['main_rgb_map', 'reg_figs',
                        'eof_summ_fig', 'raw_rgb_maps']

# write the weighted-sum EOF map to a raster file?
# NOTE: only if an extra 'save_rast' flag is fed to the call to this script
write_wt_sum_eofs_to_file = (len(sys.argv) > 2 and
                             'save_rast' in sys.argv)

# save figs?
save_it = True

# which dataset to use?
dataset = 'NIRv'

# which mask mode to use?
masking_mode = 'default'
#masking_mode = 'strict'
mask_filename_ext_dict = {'strict': '_STRICTMASK',
                          'default': ''}
mask_filename_ext = mask_filename_ext_dict[masking_mode]

# were ts standardized before EOF calculation?
standts = True
if standts:
    standts_file_substr = '_standts'
else:
    standts_file_substr = ''

# get harmonic regression design matrix
dm = phf.make_design_matrix()

# number of clusters used for global clustering?
k_global = 9


####################
# LOAD AND PREP DATA
####################

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
                                      '%s_4_EOFs_sqrt_coswts%s%s.tif') %
                         (dataset, standts_file_substr, mask_filename_ext))[:3]
eofs.rio.set_crs(4326)
# rescale each EOFs layer to [0, 1]
for i in range(eofs.shape[0]):
    eofs[i] = (eofs[i]-eofs[i].min())/(eofs[i].max()-eofs[i].min())

# load DataFrame with EOF pcts var explained and PC time series
# EOF percentages of variance explained
pc_df = pd.read_csv(os.path.join(phf.EXTERNAL_DATA_DIR,
                                 'NIRv_EOF_PC_table.csv'))
eofs_pcts = [float(c.split('_')[1].replace('pct', '')) for c in pc_df.columns]
pc_df.columns = [re.search('(?<=EOF)\d', c).group() for c in pc_df.columns]

# define focal region bounding bboxes
reg_bboxes = {
              'QLD': [1.32e7, -1.375e6, 1.373e7, -2.45e6],# far N. Queensland
              'Am': [-5.47e6, 3.1e5, -4.6e6, -0.6e6],     # mouth of Amazon
              'Ba': [-1.06e7, 4.5e6, -0.985e7, 2.87e6],   # Baja
              'GB': [-1.01e7, 5.5e6, -0.94e7, 4.6e6],     # Great Basin
              'Mad':[3.94e6, -1.4e6, 4.84e6, -3.3e6],     # Madagascar
              'Fl': [-7.525e6, 3.53e6, -7.25e6, 3.13e6],  # South Florida
              'SAf':[1.5643e6, -3.9183e6, 2.2022e6, -4.3398e6],     # S. Africa cape
              'WAu':[1.0094e7, -3.819e6, 1.1009e7, -4.374e6],  # SW Australia
              'IT': [0.55e6, 5.62e6, 1.2e6, 5.35e6],      # Po Valley, Italy
              'CB': [-8.9e6, 5.95e6, -6.95e6, 3.6e6],     # US Corn Belt
             }


# NOTE: k values determined based on manual inspection of scree plots
reg_K_vals = {
              'QLD': 3,
              'Am':  3,
              'Ba':  4,
              'GB':  3,
              'Mad': 3,
              'Fl':  3,
              'SAf': 4,
              'WAu': 4,
              'IT':  3,
              'CB':  3,
             }

# regional figs to include in main fig as keys, with values consisting of
# their labels, their gridspec indices for both the maps and line plots,
# and the positioning of their box labels within the global RGB map (expressed
# using transform=ax.transAxes)
main_fig_figsize = (20.25,13.5)
regs_in_main_gridspec_dims = (400, 600)
global_map_slices = (slice(0, 225), slice(0, 520))
regs_in_main = {'Ba': ('A.',
                       (slice(143, 373), slice(25, 125),
                        slice(375, 400), slice(30, 120)),
                       (-0.7, -1.1)),
                'GB': ('B.',
                       (slice(235, 375), slice(190, 290),
                        slice(375, 400), slice(195, 285)),
                       (-0.7, -1.1)),
                'Fl': ('C.',
                       (slice(242, 372), slice(340, 430),
                        slice(375, 400), slice(340, 430)),
                       (1.3, -1.3)),
                'Am': ('D.',
                       (slice(242, 382), slice(490, 600),
                        slice(375, 400), slice(500, 590)),
                       (1.1, -1.1)),
               }


##################
# DEFINE FUNCTIONS
##################

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


def run_clust_analysis(coeffs_rast,
                       eofs_rast,
                       reg,
                       k,
                       ax_lines,
                       k_max=12,
                       batch_size=40,
                       n_init=10,
                       max_no_improvement=10,
                       seed=1,
                       eps=1,
                       min_samples=10,
                       n_clust_neighs=50,
                      ):
    # grab all coeffs to be used for K-means clustering
    coeffs_vals = coeffs_rast.values.swapaxes(0,1).swapaxes(1,2)
    coeffs_X = coeffs_vals.reshape([np.prod(coeffs_vals.shape[:2]), 5])
    coeffs_non_nans = np.where(np.sum(np.isnan(coeffs_X), axis=1) == 0)[0]
    coeffs_X_sub = coeffs_X[coeffs_non_nans, :]
    # prep matching set of EOFs values (for color extraction later on)
    eofs_vals = eofs_rast.values.swapaxes(0,1).swapaxes(1,2)
    eofs_X = eofs_vals.reshape([np.prod(eofs_vals.shape[:2]), 3])
    eofs_X_sub = eofs_X[coeffs_non_nans, :]
    # convert to fitted phenocycle
    ts_X_sub = np.matmul(dm, coeffs_X_sub.T).T
    # standardize each phenocycle
    for i in range(ts_X_sub.shape[0]):
        ts_X_sub[i, :] = phf.standardize_array(ts_X_sub[i, :])
    # loop over K vals and make scree plot,
    # saving clustering results for the chosen k value
    # NOTE: using mini-batch K-means because it's faster
    labels = None
    centers = None
    wcss = []
    scree_k_vals = [*range(1, k_max+1)]
    for k_val in scree_k_vals:
        if seed is not None:
            np.random.seed(seed)
        labs, cents, inertia = run_kmeans_clust(ts_X_sub,
                                            k_val,
                                            batch_size,
                                            n_init,
                                            max_no_improvement,
                                           )
        # save results for the chosen k value
        if k_val == k:
            labels = deepcopy(labs)
            centers = deepcopy(cents)
        wcss.append(inertia)
    fig_scree, ax_scree = plt.subplots(1,1)
    ax_scree.plot(range(1, k_max+1), wcss)
    ax_scree.set_title(reg, fontdict={'fontsize': 24})
    fig_scree.savefig(os.path.join(phf.FIGS_DIR,
            '%s%s_LSP_clust_scree_REG_%s.png' % (dataset,
                                                  mask_filename_ext,
                                                  reg)), dpi=300)
    center_eofs = [np.median(eofs_X_sub[labels==i, :], axis=0) for i in range(k)]
    center_colors = [rgb2hex(*np.int32(ce*255)) for ce in center_eofs]
    cmap = mpl.colors.ListedColormap(center_colors)
    # create a simple map of the cluster assignments
    cluster_arr = np.ones([coeffs_X.shape[0], 1])*np.nan
    cluster_arr[coeffs_non_nans,0] = labels
    clusters = cluster_arr.reshape(coeffs_vals[:,:,0].shape)
    clusters_rxr = deepcopy(coeffs_rast[0])
    clusters_rxr[:,:] = clusters
    fig_clust_map, ax_clust_map = plt.subplots(1,1)
    phf.plot_juris_bounds(ax=ax_clust_map,
                          crs=4326,
                         )
    try:
        clusters_rxr.plot.imshow(ax=ax_clust_map,
                                 cmap=cmap, alpha=1, zorder=0)
    # ignore that annoying rioxarray error
    # ("'tuple' object has no attribute 'startswith'")
    except AttributeError as e:
        pass
    fig_clust_map.savefig(os.path.join(phf.DATA_DIR,
            '%s_%s_LSP_clust_map_REG_%s.png' % (dataset,
                                                 mask_filename_ext,
                                                 reg)), dpi=500)
    # NOTE: if this is the Great Basin, overlay the cluster raster
    #       with the cheatgrass data and run a simple ANOVA
    #       and Tukey's HSD to test differences in modeled percent annual
    #       herbaceous cover between clusters
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
        # run Monte Carlo of ANOVA (since pixel-based samples not independent)
        null_fvals = []
        null_pvals = []
        all_vals = [n for v in cluster_vals.values() for n in v]
        clust_sizes = [len(v) for v in cluster_vals.values()]
        for _ in range(1000):
            np.random.shuffle(all_vals)
            null_vals_0 = all_vals[:clust_sizes[0]]
            null_vals_1 = all_vals[clust_sizes[0]:clust_sizes[0]+clust_sizes[1]]
            null_vals_2 = all_vals[clust_sizes[0]+clust_sizes[1]:]
            assert len(null_vals_0) + len(null_vals_1) + len(null_vals_2) == len(all_vals)
            null_fval, null_pval = stats.f_oneway(null_vals_0,
                                                  null_vals_1,
                                                  null_vals_2,
                                                 )
            null_fvals.append(null_fval)
            null_pvals.append(null_pval)
        # NOTE: F-distribution is one-sided, so no need for np.abs()
        anova_permut_signif = np.mean(np.array(null_fvals)>=fval)
        # run Tukey HSD
        tuk = pairwise_tukeyhsd(
            endog=np.array([val for val_list in cluster_vals.values(
                                                ) for val in val_list]),
            groups=np.array([clust for clust, val_list in cluster_vals.items(
                                                ) for val in val_list]))
        print('\n\n==========\nGreat Basin cheatgrass statistical results:')
        print(f'\n\tANOVA: F-val: {fval}; p-val: {pval}\n\n')
        print((f"\n\tANOVA: Monte Carlo p-val: {anova_permut_signif}\n"
               f"\t(null F-vals [{np.min(null_fvals)}, {np.max(null_fvals)}])\n"
               f"\t(null P-vals [{np.min(null_pvals)}, {np.max(null_pvals)}])\n"))
        print('\tgroup means:\n\t%s' % '\t'.join(['%i: %0.2f' % (i,
                np.mean(vals)) for i, vals in cluster_vals.items()]))
        print(f"\tgroup n's:\n\t{[len(v) for v in cluster_vals.values()]}")
        print('\n\n\tTukey HSD: %s' % str(tuk))
        print(f"\n\n\tTukey HSD P-values: {tuk.pvalues}")
        print('\n\ncluster colors:\n')
        for clust, color in zip(sorted(np.unique(clusters_rxr)),
                                center_colors):
            print(f'cluster {clust}: {color}\n')
        print('\n\n==========\n\n')
    # plot cluster centers (i.e., characteristic phenocycles)
    for i in range(centers.shape[0]):
        center = centers[i, :]
        center_minmax = (center-np.min(center))/(np.max(center)-np.min(center))
        ax_lines.plot(center_minmax,
                      color=center_colors[i],
                      linewidth=3.5)
        ax_lines.set_xlim((0, 365))
        ax_lines.set_ylim(-0.055, 1.055)
    return labels, centers, center_colors


def add_phen_labs(ax, reg, text_size=13, mark_size=65, hspace_frac=0.22):
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
    elif n_lines == 5:
        hspace_top_offset = 0.25
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
    # NOTE: reordering line colors for Amazon to be first open-veg in green,
    #       then forest in blue and purple
    elif reg == 'Am':
        line_order = [line_order[1], line_order[2], line_order[0]]
    colors = [ax.get_lines()[n].get_color() for n in line_order]
    if len(colors) == 4:
        ys = np.linspace(0.98, 0.02, len(colors))
    elif len(colors) == 3:
        ys = np.linspace(0.82, 0.18, len(colors))
    for n, color in enumerate(colors):
        lab = n+1
        ax.text(385,
                ys[n]-0.08,
                #(ylim[1]-
                 #(n*hspace_frac*np.diff(ylim))-
                 #(0.5*hspace_frac)-
                 #hspace_top_offset+0),
                str(lab),
                color='black',
                size=text_size,
                clip_on=False,
               )
        ax.scatter(373,
                   ys[n],
                   #(ylim[1]-
                   # (n*hspace_frac*np.diff(ylim))-
                   # hspace_top_offset),
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


####################################
# CALCULATE WEIGHTED SUM ACROSS ITCZ
####################################

out_da_filename = ('%s_4_EOFs_sqrt_coswts%s%s'
                   '_SCALED_FOLDED_EPSG-8857.tif') % (dataset,
                                                      standts_file_substr,
                                                      mask_filename_ext,
                                                     )
out_da_filepath = os.path.join(phf.EXTERNAL_DATA_DIR, out_da_filename)

print('\n\nScaling, folding, and reprojecting EOFs...\n\n')
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
eofs_wt_sum = eofs_wt_sum.rio.write_crs(4326)
# reproject to Equal Earth
eofs_wt_sum_for_map = eofs_wt_sum.rio.reproject(8857)
# and truncate artificially inflated values along edges
# (not sure what's broken that's causing them)
for lyr in range(eofs_wt_sum_for_map.shape[0]):
    eofs_wt_sum_for_map[lyr] = eofs_wt_sum_for_map[lyr].where(
        eofs_wt_sum_for_map[lyr] < 2*eofs_wt_sum[lyr].max(), np.nan)

# mask data lying outside the edges of the
# Equal Earth projection (to prevent NZ, Sibera, etc. from wrapping around)
eofs_wt_sum_for_map = phf.mask_xarr_to_other_xarr_bbox(eofs_wt_sum_for_map,
                                                       eofs,
                                                       drop=False,
                                                       n_bbox_xcoords=4000,
                                                       n_bbox_ycoords=2000,
                                                      )
# prep ITCZ for plotting with the same map
itcz_for_plot = itcz[itcz['time'] == 'mean'].to_crs(
    eofs_wt_sum_for_map.rio.crs.to_epsg())

# if requested, write mapping-prepped, folded EOFS to file
if write_wt_sum_eofs_to_file:
    print('\n\nWriting scaled, folded, reprojected EOFs to file...\n\n')
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
    out_da.rio.to_raster(out_da_filepath)


#######################
## PLOT EOF SUMMARY FIG
#######################

if what_to_plot == 'eof_summ_fig':
    # create EOF fig
    fig_eof = plt.figure(figsize=(14,20))
    gs = GridSpec(4, 3, figure=fig_eof, height_ratios=[1,1,1,1.1])
    # make maps of EOFS 1 through 4, all raw
    # NOTE: reloading, to grab all 4
    eofs = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                          '%s_4_EOFs_sqrt_coswts%s%s.tif') %
                        (dataset, standts_file_substr, mask_filename_ext))[:4]
    # rescale each EOFs layer to [0, 1]
    for i in range(eofs.shape[0]):
        eofs[i] = (eofs[i]-eofs[i].min())/(eofs[i].max()-eofs[i].min())
    eofs.rio.set_crs(4326)
    eofs_for_map = eofs.rio.reproject(8857)
    # mask to Equal Earth projection bounds
    eofs_for_map = phf.mask_xarr_to_other_xarr_bbox(eofs_for_map,
                                                    eofs,
                                                    drop=False,
                                                    n_bbox_xcoords=4000,
                                                    n_bbox_ycoords=2000,
                                                   )
    for i in range(4):
        ax_map = fig_eof.add_subplot(gs[i, :2])
        ax_pc = fig_eof.add_subplot(gs[i, 2])
        # add colorbar axes at bottom
        if i == 3:
            add_colorbar=True
            divider = make_axes_locatable(ax_map)
            cbar_ax = divider.append_axes('bottom', size='7%', pad=0.2)
            cbar_ax.tick_params(labelsize=12)
            cbar_kwargs = {'orientation': 'horizontal'}
        else:
            add_colorbar=False
            cbar_ax = None
            cbar_kwargs = None
        eofs_for_map[i] = eofs_for_map[i].where(
                eofs_for_map[i] < 2*eofs[i].max(), np.nan)
        # reset the 'long_name' attr, to avoid the stupid "AttributeError:
        # 'tuple' object has no attribute 'startswith'" that is thrown by
        # rioxarry if I otherwise try to use add_colorbar=True
        eofs_for_map.attrs['long_name'] = ''
        eofs_for_map[i].plot.imshow(ax=ax_map,
                                    cmap='coolwarm',
                                    add_colorbar=add_colorbar,
                                    cbar_ax=cbar_ax,
                                    cbar_kwargs=cbar_kwargs,
                                    alpha=1,
                                    zorder=0,
                                   )
        if i == 3:
            cbar_ax.set_xlabel('EOF', fontdict={'size': 20})
        phf.plot_juris_bounds(ax_map,
                              lev0_linewidth=0.2,
                              lev0_alpha=0.7,
                              lev1_linewidth=0.1,
                              lev1_alpha=0.5,
                              strip_axes=True,
                             )
        ax_map.set_ylim(eofs_for_map.rio.bounds()[1::2])
        ax_map.text(0.92*ax_map.get_xlim()[0], 0.92*ax_map.get_ylim()[0],
                    'EOF %i:\n%0.2f%%' % (i+1, eofs_pcts[i]),
                    fontdict={'fontsize': 26})
        ax_map.set_aspect('equal')
        phf.set_upper_ylim(ax_map)
        # plot PC time series
        ax_pc.plot(pc_df.loc[:, str(i)], linewidth=2, color='black')
        ax_pc.set_xlabel('day of year', fontdict={'fontsize': 20})
        ax_pc.set_ylabel('PC', fontdict={'fontsize': 20})
        ax_pc.tick_params(axis='both', labelsize=12)
        # square the PC axes
        ax_pc.set_box_aspect(1)
    del eofs_for_map
    fig_eof.subplots_adjust(left=0.02, right=0.98, bottom=0, top=1, hspace=0)
    fig_eof.subplots_adjust(wspace=.3)
    if save_it:
        fig_eof.savefig(os.path.join(phf.FIGS_DIR,
                                     'FIG_SUPP_%s_EOF_summary%s.pdf' % (dataset,
                                                mask_filename_ext)), dpi=600)


#############################
# PLOT UNTRANSFORMED RGB MAPS
#############################

if what_to_plot == 'raw_rgb_maps':
    # create untransformed RGB figure
    eofs_for_map = eofs.rio.write_crs(4326)
    eofs_for_map = eofs_for_map.rio.reproject(8857)
    # mask to bounds of Equal Earth projection
    # mask to Equal Earth projection bounds
    eofs_for_map = phf.mask_xarr_to_other_xarr_bbox(eofs_for_map,
                                                    eofs,
                                                    drop=False,
                                                    n_bbox_xcoords=4000,
                                                    n_bbox_ycoords=2000,
                                                   )
    fig_untrans = plt.figure(figsize=(20,20))
    ax_top = fig_untrans.add_subplot(2,1,1)
    ax_bot = fig_untrans.add_subplot(2,1,2)
    axs = {0: ax_top, 1: ax_bot}
    for row in range(2):
        ax = axs[row]
        if row == 0:
            eofs_for_map.plot.imshow(ax=ax, zorder=0)
        elif row == 1:
            eofs_trans = deepcopy(eofs_for_map)
            # NOTE: folding all 3 EOFs, on the assumption that a hemispheric signal is
            #       embedded in each, rather than eyeballing whether or not one is...
            for i in range(3):
                eofs_trans[i,:,:] = 1 - eofs_trans[i,:,:]
            eofs_trans = eofs_trans.where(np.abs(eofs_trans)<1e37)
            eofs_trans.plot.imshow(ax=ax, zorder=0)
        phf.plot_juris_bounds(ax,
                              lev0_linewidth=0.2,
                              lev0_alpha=0.7,
                              lev1_linewidth=0.1,
                              lev1_alpha=0.5,
                              strip_axes=True,
                              reset_axlims=False,
                             )
        itcz_for_plot.plot(ax=ax,
                      linestyle='--',
                      linewidth=0.5,
                      color='k',
                      alpha=0.5,
                     )
        ax.set_ylim(eofs_for_map.rio.bounds()[1::2])
        ax.set_aspect('equal')
        phf.set_upper_ylim(ax)
    fig_untrans.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98,
                                hspace=0.05)
    if save_it:
        fig_untrans.savefig(os.path.join(phf.FIGS_DIR,
                                         'FIG_SUPP_untransformed_EOF_maps.pdf'), dpi=600)

    del eofs_for_map



###################
# PLOT MAIN RGB MAP
###################

if what_to_plot == 'main_rgb_map':

    # create main figure
    fig_1 = plt.figure(figsize=main_fig_figsize)
    gs = GridSpec(*regs_in_main_gridspec_dims, figure=fig_1)
    ax_map = fig_1.add_subplot(gs[global_map_slices[0], global_map_slices[1]])

    # get rid of first 198 columns and first 310 rows, to reduce white space
    # (first column with unmasked values is 199, and row is 312)
    eofs_wt_sum_for_map_trunc = eofs_wt_sum_for_map[:, 310:,198:]
    # plot the RGB map of EOFs 1-3 together
    eofs_wt_sum_for_map_trunc.plot.imshow(ax=ax_map,
                                          add_colorbar=False,
                                          alpha=1,
                                          zorder=0,
                                         )
    phf.plot_juris_bounds(ax_map,
                          lev0_linewidth=0.2,
                          lev0_alpha=0.7,
                          lev1_linewidth=0.1,
                          lev1_alpha=0.5,
                          strip_axes=True,
                         )
    # add the ITCZ as a thin gray line
    itcz_for_plot.plot(ax=ax_map,
                      linestyle='--',
                      linewidth=0.5,
                      color='k',
                      alpha=0.5,
                     )
    ax_map.set_ylim(eofs_wt_sum_for_map_trunc.rio.bounds()[1::2])
    ax_map.set_aspect('equal')
    phf.set_upper_ylim(ax_map)
    ax_map.spines['bottom'].set_color('none')
    ax_map.spines['top'].set_color('none')
    ax_map.spines['right'].set_color('none')
    ax_map.spines['left'].set_color('none')

    # read in and plot the global ts clusters' time series
    global_clust_cents_path = os.path.join(phf.EXTERNAL_DATA_DIR,
                                           f'NIRv_LSP_{k_global}_clusts_cents_GLOBAL.pkl')
    with open(global_clust_cents_path, 'rb') as f:
        clust_centers = pickle.load(f)
    # read in the cluster-assignment array and use that to calculate each
    # cluster's color as its median value in the ITCZ-folded EOFs raster
    global_clust_arr_path = os.path.join(phf.EXTERNAL_DATA_DIR,
                                         f'NIRv_LSP_{k_global}_clusts_arr_GLOBAL.txt')
    global_clust_arr = np.loadtxt(global_clust_arr_path)
    clust_colors = []
    for k in range(k_global):
        eofs_masked = eofs_wt_sum.where(global_clust_arr==k)
        color = np.nanmedian(eofs_masked.values.reshape(3, -1), axis=1)
        clust_colors.append(color)
    # order by doy of max value
    order_top_bot = np.argsort([np.argmax(cent) for cent in clust_centers])
    otb = order_top_bot
    order_top_bot = otb[[2, 4, 3, 5, 8, 7, 6, 1, 0]]
    for n, i in enumerate(order_top_bot):
        cent = clust_centers[i]
        color = clust_colors[i]
        ax_lines_slice_x = slice(510, 595)
        ax_lines_slice_y = slice(n*25, (n+1)*25)
        ax_lines = fig_1.add_subplot(gs[ax_lines_slice_y,
                                        ax_lines_slice_x])
        ax_lines.plot(cent, color=color, linewidth=3.5)
        phf.strip_axes_labels_and_ticks(ax_lines)
        ax_lines.set_xlabel('')
        ax_lines.set_ylabel(f'  {n+1}', size=16, rotation=0)
        ax_lines.yaxis.set_label_position("right")
        ax_lines.set_xlim(0, 364)
        if n == k_global-1:
            # TODO: ADD SOUTHERN HEM LABELS AS SECOND SET OF LABELS ROTATED 180?
            ax_lines.set_xticks(np.linspace(0, 365, 5),
                          ['Jan\n—\nJul',
                           'Apr\n—\nOct',
                           'Jul\n—\nJan',
                           'Oct\n—\nApr',
                           'Jan\n—\nJul',
                          ])
        # adjust focal-region plotting and add bounding boxes
        for month in np.linspace(0, 365, 13):
            ax_lines.axvline(month, color='black',
                       linewidth=0.25, linestyle=':', alpha=0.75, zorder=0)
        ax_lines.tick_params(labelsize=13, rotation=0)
        for axis in ['top','bottom','left','right']:
            ax_lines.spines[axis].set_linewidth(2)
            ax_lines.spines[axis].set_color('black')

    # also create identical miniature figure, for embedding in top left of
    # supplemental figs
    asp_ratio = eofs_wt_sum_for_map_trunc.shape[1]/eofs_wt_sum_for_map_trunc.shape[2]
    fig_miniglobal = plt.figure(figsize=(6,6*asp_ratio))
    ax_mini = fig_miniglobal.add_subplot(1,1,1)
    eofs_wt_sum_for_map_trunc.plot.imshow(ax=ax_mini,
                                          add_colorbar=False,
                                          alpha=1,
                                          zorder=0,
                                         )
    phf.plot_juris_bounds(ax_mini,
                          lev0_linewidth=0.2,
                          lev0_alpha=0.7,
                          lev1_linewidth=0.1,
                          lev1_alpha=0.5,
                          strip_axes=True,
                         )
    # add the ITCZ as a thin gray line
    itcz_for_plot.plot(ax=ax_mini,
                      linestyle='--',
                      linewidth=0.5,
                      color='k',
                      alpha=0.5,
                     )
    ax_mini.set_ylim(eofs_wt_sum_for_map_trunc.rio.bounds()[1::2])
    phf.set_upper_ylim(ax_mini)
    ax_mini.spines['bottom'].set_color('none')
    ax_mini.spines['top'].set_color('none')
    ax_mini.spines['right'].set_color('none')
    ax_mini.spines['left'].set_color('none')
    ax_mini.set_aspect('equal')
    fig_miniglobal.subplots_adjust(left=0.0,
                                   right=0.99,
                                   bottom=0.04,
                                   top=0.98)
    if save_it:
        fig_miniglobal.savefig(os.path.join(phf.FIGS_DIR, 'SUBFIG_mini_global_LSP_RGB_map.png'),
                      dpi=400)



########################
# PLOT REGIONAL RGB MAPS
########################

if what_to_plot in ['reg_figs', 'main_rgb_map']:

    ####################
    # PLOT FOCAL REGIONS
    #################### 

    for reg, bbox in reg_bboxes.items():
        # plot in main fig?
        in_main = reg in regs_in_main
        if ((in_main and what_to_plot == 'main_rgb_map') or
            (not in_main and what_to_plot == 'reg_figs')):
            print('now plotting region: %s' % reg)
            if not in_main:
                # create regional figure
                dims = get_fig_dims_from_reg_bbox(bbox)
                fig_reg, ax_reg = plt.subplots(1, figsize=dims)
                fig_lines, ax_lines = plt.subplots(1, figsize=(4,1))
            else:
                label, gs_slices, box_lab_coords = regs_in_main[reg]
                ax_reg = fig_1.add_subplot(gs[gs_slices[0], gs_slices[1]])
                ax_lines = fig_1.add_subplot(gs[gs_slices[2], gs_slices[3]])
            # plot the focal region
            eofs_wt_sum_for_map.sel(x=slice(bbox[0], bbox[2]),
                                    y=slice(bbox[1], bbox[3])).plot.imshow(ax=ax_reg,
                                                                           zorder=0,
                                                                    add_colorbar=False,
                                                                          )
            phf.plot_juris_bounds(ax_reg, lev0_alpha=0.8, strip_axes=True)
            ax_reg.set_xlim(bbox[0], bbox[2])
            ax_reg.set_ylim(bbox[3], bbox[1])
            ax_reg.set_aspect('equal')
            # subset the focal region's data and run clustering
            eofs_wt_sum_foc = eofs_wt_sum_for_map.sel(x=slice(bbox[0], bbox[2]),
                                                      y=slice(bbox[1], bbox[3]),
                                                     ).rio.write_crs(8857)
            eofs_wt_sum_foc = eofs_wt_sum_foc.rio.reproject(4326)
            eofs_wt_sum_foc = eofs_wt_sum_foc.where(eofs_wt_sum_foc !=
                                                    eofs_wt_sum_foc._FillValue, np.nan)
            coeffs_foc = coeffs.rio.reproject(8857).sel(x=slice(bbox[0], bbox[2]),
                                                        y=slice(bbox[1], bbox[3]),
                                                     ).rio.write_crs(8857)
            coeffs_foc = coeffs_foc.rio.reproject(4326)
            coeffs_foc = coeffs_foc.where(coeffs_foc != coeffs_foc._FillValue, np.nan)
            # run K-means clustering
            K = reg_K_vals[reg]
            (clust_labs,
             clust_cents,
             clust_colors) = run_clust_analysis(coeffs_rast=coeffs_foc,
                                                eofs_rast=eofs_wt_sum_foc,
                                                reg=reg,
                                                k=K,
                                                ax_lines=ax_lines,
                                               )
            # add line plot labels
            add_phen_labs(ax_lines, reg)
            # adjust focal-region plotting and add bounding boxes
            phf.strip_axes_labels_and_ticks(ax_lines)
            ax_lines.set_xticks(np.linspace(0, 365, 5))
            ax_lines.set_xticklabels(['Jan', 'Apr', 'Jul', 'Oct', 'Jan'])
            for month in np.linspace(0, 365, 13):
                ax_lines.axvline(month, color='black',
                                 linewidth=0.25, linestyle=':', alpha=0.75, zorder=0)
            ax_lines.tick_params(labelsize=13, rotation=0)
            for axis in ['top','bottom','left','right']:
                ax_lines.spines[axis].set_linewidth(2)
            if in_main:
                # black boxes around regional subplots within main fig
                for axis in ['top','bottom','left','right']:
                    ax_reg.spines[axis].set_color('black')
                    ax_lines.spines[axis].set_color('black')
                # label subplots, and add matching-label boxes to the global map
                box_for_map = box(*bbox)
                gdf_for_map = gpd.GeoDataFrame(pd.DataFrame({'i': [0],
                                                             'geometry': box_for_map}),
                                               crs=8857,
                                               geometry='geometry',
                                              )
                gdf_for_map.plot(ax=ax_map,
                                 color='none',
                                 edgecolor='black',
                                 linewidth=1.5,
                                )
                ax_map.text(x=bbox[0] + box_lab_coords[0]*(bbox[2]-bbox[0]),
                            y=bbox[3] + box_lab_coords[1]*(bbox[3]-bbox[1]),
                            s=label,
                            size=14,
                            weight='bold',
                           )
                ax_reg.text(x=-0.2,
                            y=1.0,
                            s=label,
                            size=30,
                            weight='bold',
                            clip_on=False,
                            transform=ax_reg.transAxes,
                           )

            # if this is for a supplemental regional fig then save components separately
            if not in_main:
                fig_reg.subplots_adjust(left=0.01,
                                        right=0.99,
                                        bottom=0.04,
                                        top=0.98,
                                       )
                fig_lines.subplots_adjust(left=0.04,
                                          right=0.94,
                                          bottom=0.24,
                                          top=0.97,
                                         )
                if save_it:
                    fig_reg.savefig(os.path.join(phf.FIGS_DIR, 'SUBFIG_REG_%s_%s_RGB_EOF_reg_map%s.png' % (reg, dataset,
                                        mask_filename_ext)), dpi=700)
                    fig_lines.savefig(os.path.join(phf.FIGS_DIR, 'SUBFIG_REG_%s_LINES_%s_RGB_EOF_reg_map%s.png' % (reg, dataset,
                                        mask_filename_ext)), dpi=700)

# save main figure after regional maps have been added
if what_to_plot == 'main_rgb_map':
    # save figure
    fig_1.subplots_adjust(left=0.0,
                          right=0.99,
                          bottom=0.04,
                          top=0.98)
    if save_it:
        fig_1.savefig(os.path.join(phf.FIGS_DIR, 'FIG_1_LSP_RGB_map.pdf'),
                      dpi=700)



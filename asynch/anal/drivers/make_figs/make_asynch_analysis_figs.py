import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from matplotlib.patches import Polygon as mplPolygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as path_effects
from collections import Counter as C
from shapely.geometry import Polygon as shapelyPolygon
from math import pi
from copy import deepcopy
import seaborn as sns
import statsmodels.api as sm
import palettable
import cmocean
import colorsys
import cmocean
import dask
import sys
import re
import os

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf


# which plots to create?
what_to_plot = sys.argv[1].strip().lower()
assert what_to_plot in ['main', 'error_supp', 'predom_supp']
if what_to_plot == 'main':
    only_top_covars = True
else:
    only_top_covars = False


# set data-directory path
data_dir = phf.EXTERNAL_RF_DATA_DIR

# general plotting params:
title_fontsize = 12
rowlab_fontsize = 18
axislab_fontsize = 13
ticklab_fontsize = 12
cbarlab_fontsize = 12
cbar_axlab_fontsize = 35
cbar_ticklab_fontsize = 24
annot_fontsize = 14
scat_label_fontsize = 10
scat_label_fontcolor = 'black'
scat_label_linewidth = 3
fig_width = 10.5
fig_height = 12
dpi = 400
n_ticklabels = 5
subplots_adj_left=0.06
subplots_adj_bottom=0.05
subplots_adj_right=0.96
subplots_adj_top=0.95
subplots_adj_wspace=0.14
subplots_adj_hspace=0.30
central_curve_color = '#220000'
central_linewidth = 2
central_alpha = 1
neighbor_linewidth = 0.1
neighbor_alpha = 0.5
rad_linestyle = ':'
rad_color = 'white'
rad_linewidth = 2
rad_alpha = 0.75
rad_mask_alpha = 0.175
supp_axlabel_fontdict = {'fontsize': 18}
cbarlabel_fontdict = {'fontsize': 14}
cbarticklabel_fontsize = 12
partlabel_fontsize=40
plot_crs = 8857


# asynch params for main results:
# main variable
var = 'NIRv'
# asynch neighborhood radius
neigh_rad = 100
# whether or not to use results of the model that included geocoods as covars
include_coords = 'y'

# asynch percentile threshold below which to drop pixels from summary map
asynch_thresh = 75

# process the data raw, or read in processed files?
if what_to_plot == 'main':
    process_raw = not os.path.isfile(os.path.join(data_dir,
                                                  'SHAP_predom_top.tif'))
else:
    process_raw = not os.path.isfile(os.path.join(data_dir,
                                                  'SHAP_predom_all.tif'))

# asynch colormap
#asynch_cmap = 'viridis'
asynch_cmap = palettable.scientific.sequential.LaJolla_20_r.mpl_colormap

# and predominance colormap
if only_top_covars:
    #predom_cmap = palettable.scientific.diverging.Lisbon_20.get_mpl_colormap()
    predom_cmap = palettable.scientific.diverging.Berlin_20.get_mpl_colormap()
else:
    # NOTE: adapted from 'light' colormap at https://personal.sron.nl/~pault/#sec:qualitative
    colors = [
              '#77AADD', # light blue -> ppt.asy
              '#FFAABB', # pink -> tmp.min.asy
              '#BBCC33', # pear -> veg.ent
              '#EE8866', # orange -> tmp.max.asy
              '#EEDD88', # light yellow -> cld.asy
              '#99DDFF', # light cyan -> def.asy
              '#44BB99', # mint -> vrm.med
              #'#AAAA00', # olive
              #'#ADADAD', # light grey
             ]
    predom_cmap = ListedColormap(colors)


# set up SHAP data:
def minmax_scale(vals, bot_pctile, top_pctile):
    """
    min-max scaler
    """
    scaled = (vals - np.nanpercentile(vals, bot_pctile))/(
                     np.nanpercentile(vals, top_pctile) -
                     np.nanpercentile(vals, bot_pctile))
    # clip to [0,1]
    scaled = np.clip(scaled, 0, 1)
    return scaled


# top covariates (in case only mapping those)
top_covars = ['ppt.asy', 'tmp.min.asy']

# NOTE: THIS DICT AND THE FOLLOWING ARE ORDERED ACCORDING TO DECREASING
#       SHAP-BASED COVAR IMPORTANCE IN OUR MAIN RF MODEL
covar_cbar_labels = {
    'ppt.asy': '$\Delta dist_{seas_{P}}/\Delta  dist_{geo}$',
    'tmp.min.asy': '$\Delta dist_{seas_{T_{min}}}/\Delta  dist_{geo}$',
    'veg.ent': '$entropy$',
    'tmp.max.asy': '$\Delta dist_{seas_{T_{max}}}/\Delta  dist_{geo}$',
    'cld.asy': '$\Delta dist_{seas_{cld}}/\Delta  dist_{geo}$',
    'def.asy': '$\Delta dist_{seas_{def}}/\Delta  dist_{geo}$',
    'vrm.med': '$med_{VRM}$',
}

covar_longnames = {
    'ppt.asy': 'precipitation\nasynchrony',
    'tmp.min.asy': 'min. temperature\nasynchrony',
    'veg.ent': 'vegetation\nentropy',
    'tmp.max.asy': 'max. temperature\nasynchrony',
    'cld.asy': 'cloud\nasynchrony',
    'def.asy': 'CWD\nasynchrony',
    'vrm.med': 'median vector\nruggedness metric',
    'codom':   'codominance',
}

cbar_axlab_dict = {'NIRv main': 'LSP asynchrony',
                   'NIRv': '$NIR_{V}\ asynch\ (\Delta NIR_{V}/\Delta m)$',
                   'SIF': '$SIF\  asynch\ (\Delta (mW m^{-2} sr^{-1} nm^{-1})/\Delta m)$',
                   'tmmn': '$tmp_{min}\ asynch\ (\Delta ^{\circ} C/\Delta m)$',
                   'tmmx': '$tmp_{max}\ asynch\ (\Delta ^{\circ} C/\Delta m)$',
                   'pr': '$ppt\ asynch\ (\Delta mm/\Delta m)$',
                   'def': '$cwd\ asynch\ (\Delta mm/\Delta m)$',
                   'cloud': '$cloud\ asynch\ (\Delta \%\ cover/\Delta m)$',

                  }

def plot_juris_bounds(ax, black_zorder=0, subnat_zorder=1, nat_zorder=2,
                      polys_color='#060606',
                     ):
    """
    plot national and subnational jurisdictional bounds
    """
    if black_zorder is not None:
        phf.plot_juris_bounds(ax,
                              lev1=False,
                              lev0_color=polys_color,
                              lev0_linewidth=0,
                              lev0_alpha=0.2,
                              lev0_zorder=black_zorder,
                              crs=plot_crs,
                              strip_axes=False,
                             )
    phf.plot_juris_bounds(ax,
                          lev1_linecolor='gray',
                          lev1_linewidth=0.3,
                          lev1_alpha=0.7,
                          lev1_zorder=subnat_zorder,
                          lev0_linecolor='gray',
                          lev0_linewidth=0.5,
                          lev0_alpha=0.8,
                          lev0_zorder=nat_zorder,
                          crs=plot_crs,
                          strip_axes=False,
                         )


def map_asynch(fig, cbar_axlab,
               gs=None, main_fig=True, var='NIRv',
               cbar_axlab_fontsize=cbar_axlab_fontsize,
               cbar_ticklab_fontsize=cbar_ticklab_fontsize):
    """
    plot the a phenological or climatic asynchrony map
    """
    assert var in ['NIRv', 'SIF', 'tmmn', 'tmmx', 'pr', 'def', 'cloud']
    if var in ['NIRv', 'SIF']:
        files = [f for f in os.listdir(phf.EXTERNAL_DATA_DIR) if
                                        re.search('%s_STRICT_asynch' % var, f)]
    else:
        files = [f for f in os.listdir(phf.EXTERNAL_DATA_DIR) if
                                        re.search('%s_asynch' % var, f)]
    # cut down to just one file, if this is for the main fig
    if main_fig:
        files = [f for f in files if re.search('asynch_100km', f)]
    # otherwise arrange in top-down order of increasing neighborhood radius
    else:
        reordered_files = []
        for neigh_rad in [50, 100, 150]:
            neigh_rad_file = [f for f in files if re.search('asynch_%ikm' %
                                                            neigh_rad, f)]
            assert len(neigh_rad_file) == 1
            reordered_files.append(neigh_rad_file[0])
        files = reordered_files
    for ax_ct, file in enumerate(files):
        # get the neighborhood radius
        neigh_rad = int(re.search('(?<=asynch_)\d{2,3}(?=km\.tif)',
                                  file).group())
        # either grab the upper half of the main fig
        if main_fig:
            ax = fig.add_subplot(gs[:130, :])
        # or grab the next row of the supp fig
        else:
            ax = fig.add_subplot(3,1,ax_ct+1)
        # partition off a separate axis for the colormap
        divider = make_axes_locatable(ax)
        if main_fig:
            where = 'bottom'
            orientation = 'horizontal'
            size = '7%'
        else:
            where = 'right'
            orientation = 'vertical'
            size = '4%'
        cax = divider.append_axes(where, size=size, pad=0.2)
        plot_juris_bounds(ax, 0, 2, 3)
        # read in the raster data and prepare it
        rast = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                              file), masked=True)[0]
        rast_proj = rast.rio.write_crs(4326).rio.reproject(plot_crs)
        # mask outside global Equal Earth bounds
        rast_proj = phf.mask_xarr_to_other_xarr_bbox(rast_proj, rast)
                # NOTE: annoying AttributeError is because da.attrs['long_name']
        #       is retained as a tuple of names (rather than being subsetted
        #       by indexing) when I index a single layer out of an
        #       xarray.core.dataarray.DataArray;
        #       for now, a hacky fix is just assigning a string to that attr
        rast_proj.attrs['long_name'] = ''
        rast_proj.plot.imshow(ax=ax,
                         cmap=asynch_cmap,
                         vmin=np.nanpercentile(rast_proj, 1),
                         vmax=np.nanpercentile(rast_proj, 99),
                         add_colorbar=True,
                         cbar_ax=cax,
                         cbar_kwargs = {'orientation': orientation},
                         zorder=1,
                        )
        cax.tick_params(labelsize=cbar_ticklab_fontsize)
        if main_fig:
            cax.set_xlabel(cbar_axlab, fontdict={'fontsize': cbar_axlab_fontsize})
        else:
            cax.set_ylabel(cbar_axlab, fontdict={'fontsize': cbar_axlab_fontsize})
                # format axes
        ax.set_xlim(rast_proj.rio.bounds()[0::2])
        ax.set_ylim(rast_proj.rio.bounds()[1::2])
        # NOTE: chopping off western edge because the equal earth projection
        #       makes NZ appear twice
        ax.set_xlim(0.95 * ax.get_xlim()[0], ax.get_xlim()[1])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks(())
        ax.set_yticks(())
        # add axis title, if not main figure
        if main_fig:
            ax.set_title('')
        else:
            ax.set_title('%i km neighborhood' % neigh_rad,
                         fontdict={'fontsize': 21})
        del rast_proj, rast


def map_predom(ax, predom,
               xlim, ylim,
               covars, which_covars='top',
               cbar_axlab='covariate',
               cbar_axlab_fontsize=cbar_axlab_fontsize,
               cbar_ticklab_fontsize=cbar_ticklab_fontsize,
               polys_color='#060606',
              ):
    """
    plot a summary map of RF covariate predominance
    """
    divider = make_axes_locatable(ax)
    if which_covars == 'top':
        cax = divider.append_axes('bottom', size='7%', pad=0.2)
        cbar_kwargs = {'orientation': 'horizontal'}
        vmin = -0.99
        vmax = 0.99
    else:
        cax = cbar_kwargs = vmin = vmax = None
    predom.plot.imshow(ax=ax,
                       cmap=predom_cmap,
                       add_colorbar=(which_covars == 'top'),
                       cbar_ax=cax,
                       cbar_kwargs=cbar_kwargs,
                       zorder=1,
                       vmin=vmin,
                       vmax=vmax,
                      )
    ax.set_xticks(())
    ax.set_yticks(())
    plot_juris_bounds(ax, 0, 5, 6, polys_color=polys_color)
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('')
    # set axes limits correctly
    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        ax.set_ylim(predom.rio.bounds()[1::2])
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_xlim(predom.rio.bounds()[0::2])
        # NOTE: chopping off western edge because the equal earth projection
        #       makes NZ appear twice
        ax.set_xlim(0.95 * ax.get_xlim()[0], ax.get_xlim()[1])
    # make manual colorbar
    if which_covars == 'all':
        polys = []
        tick_locs = []
        tick_labs = []
        cax = divider.append_axes('bottom', size='7%', pad=0.2)
        for i, covar in enumerate(covars):
            polys.append(mplPolygon([(i, 0),
                                  (i, 1),
                                  (i+1, 1),
                                  (i+1, 0),
                                  (i, 0)]))
            tick_locs.append(i+0.5)
            tick_labs.append(covar_longnames[covar])
        patches = PatchCollection(polys, alpha=1, edgecolor='k')
        patches.set_color([mpl.colors.hex2color(c) for c in colors])
        cax.add_collection(patches)
        cax.set_xlim((0, len(covars)))
    else:
        cbar_covars = [covars[0], 'codom', covars[1]]
        tick_locs = [-0.9, 0, 0.9]
        tick_labs = [covar_longnames[covar] for covar in cbar_covars]
    cax.set_ylim((0,1))
    cax.set_xlabel(cbar_axlab, fontdict={'fontsize': cbar_axlab_fontsize})
    cax.set_xticks(tick_locs)
    cax.set_xticklabels(tick_labs, fontdict={'fontsize':
                                                 cbar_ticklab_fontsize})
    cax.set_ylabel('')
    cax.set_yticks(())
    if which_covars == 'top':
        cax.tick_params(axis=u'both', which=u'both',length=0)



if __name__ == '__main__':

    if what_to_plot in ['main', 'predom_supp']:

    ###################
    # CREATE PREDOM FIG
    ###################

        # stack all covars' SHAP maps
        files = [f for f in os.listdir(data_dir) if re.search(
            f'SHAP_map_{include_coords}COORDS_[cdptv].*_{var}_{neigh_rad}km.tif', f)]
        # reorder per decreasing SHAP covar importance in our main RF model
        # (i.e., the order of the covar_cbar_labels dict)
        reordered_files = []
        for covar in covar_cbar_labels.keys():
            for f in files:
                if re.search(f'{covar}_{var}', f):
                    reordered_files.append(f)
        all_files = reordered_files
        shap_maps = []
        covars = []
        if only_top_covars:
            files = [f for f in all_files if re.search(f'(?<=COORDS_).*(?=_{var})', f).group() in top_covars]
            for f in files:
                covar = re.search(f'(?<=COORDS_).*(?=_{var})', f).group()
                covars.append(covar)
                shap_map = rxr.open_rasterio(os.path.join(data_dir, f), masked=True)[0]
                shap_map.name = covar
                shap_maps.append(shap_map)
            xr_stack = xr.concat(shap_maps, dim=pd.Index(covars, name='covar'))
            da = dask.array.stack(xr_stack)

        else:
            files = all_files
            for f in files:
                covar = re.search(f'(?<=COORDS_).*(?=_{var})', f).group()
                covars.append(covar)
                shap_map = rxr.open_rasterio(os.path.join(data_dir, f), masked=True)[0]
                shap_map.name = covar
                shap_maps.append(shap_map)
            xr_stack = xr.concat(shap_maps, dim=pd.Index(covars, name='covar'))
            da = dask.array.stack(xr_stack)

        # rechunk (to run slowly but successfully on laptop)
        # NOTE: axis sizes apparently needn't be exact multiples of their chunk sizes!
        da = da.rechunk(chunks=(da.shape[0], 53, 500))

        if process_raw:
            print('\n\nPROCESSING RAW DATA; PLEASE HOLD...\n\n')
            if only_top_covars:
                # for main figure, calculate a normalized-difference predominance index
                # of the absolute SHAP values
                # NOTE: WHY?
                    # 1. we only want to know about only the relative magnitudes of the
                    #    modeled effects of the two top covariates, not their directional
                    #    effects
                    # 2. we want to normalize each pixel to its own range of values,
                    #    such that predominance is equivalent at pixels with values
                    #    (for the 2 covars) of [0, 0.5], [0, 1], [0, 2.5]
                    #    as well as at pixels with values of [0.1, 0.1], [0.5, 0.5],
                    #    [1.5, 1.5], etc...
                # NOTE: subtracting the top covar (asy.ppt) from the second covar
                #       (asy.tmp.min) gives us a result for which negative values
                #       indicate asy.ppt predominance and positive values asy.tmp.min
                #       predominance (and thus blue and pink colors respectively, which
                #       will match the first two colors used in the all-covars map)
                predom = deepcopy(shap_maps[0])
                predom.values = ((np.abs(da[1]) - np.abs(da[0]))/
                                 (np.abs(da[1]) + np.abs(da[0]))).compute()
                file_suffix = 'top'
            else:
                # hue: predominant var
                predom = deepcopy(shap_maps[0])
                predom.values = np.argmax(np.abs(da), axis=0).compute()
                file_suffix = 'all'

            # save to disk
            predom.rio.to_raster(os.path.join(data_dir,
                                              f'SHAP_predom_{file_suffix}.tif'))


        else:
            if only_top_covars:
                file_suffix = 'top'
            else:
                file_suffix = 'all'
            predom = rxr.open_rasterio(os.path.join(data_dir,
                            f'SHAP_predom_{file_suffix}.tif'), masked=True)[0]

            # reproject rasters
            predom_proj = predom.rio.write_crs(4326).rio.reproject(plot_crs, nodata=np.nan)
            predom = phf.mask_xarr_to_other_xarr_bbox(predom_proj, predom)

            # asynchrony will be used to mask out lower-asynchrony pixels
            asynch = minmax_scale(rxr.open_rasterio(os.path.join(data_dir,
                    f'{var}_STRICT_asynch_{neigh_rad}km.tif'), masked=True)[0], 0, 99)
            # mask out pixels less than the percentile threshold
            asynch = asynch.where(asynch >= np.nanpercentile(asynch, asynch_thresh))
            # reproject and rescale
            asynch_reproj = asynch.rio.reproject_match(predom)
            asynch = phf.mask_xarr_to_other_xarr_bbox(asynch_reproj, asynch)
            asynch = minmax_scale(asynch, 1, 99)
            assert np.nanmin(asynch) == 0 and np.nanmax(asynch) == 1

            # mask other two layers to the thresholded asynch map
            predom = predom.where(pd.notnull(asynch))



        if what_to_plot == 'main':
            # set up main figure
            fig_main = plt.figure(figsize=(24, 26))
            gs = fig_main.add_gridspec(ncols=270, nrows=270)


            # plot the main asynch map
            map_asynch(fig_main, cbar_axlab_dict['NIRv main'], gs=gs, main_fig=True, var='NIRv')

            # plot the RF summary map
            ax = fig_main.add_subplot(gs[140:,:])
            map_predom(ax, predom,
                       xlim=fig_main.axes[0].get_xlim(),
                       ylim= fig_main.axes[0].get_ylim(),
                       covars=covars,
                       which_covars='top',
                       cbar_axlab='predominant driver of LSP asynchrony',
                      )

            # finish plot formatting
            fig_main.subplots_adjust(left=0.05,
                                     right = 0.95,
                                     bottom=0.03,
                                     top=0.98,
                                     hspace=0.8,
                                     wspace=0.0
                                    )

            # add labels for parts A. and B.
            for i, lett in zip([0, 2], ['A.', 'B.']):
                fig_main.axes[i].text(fig_main.axes[i].get_xlim()[0]*0.98,
                                      fig_main.axes[i].get_ylim()[1]*0.8,
                                      lett,
                                      size=60,
                                      weight='bold',
                                     )

            fig_main.savefig('FIG_asynch_and_rf_results.png', dpi=500)

            del fig_main

        elif what_to_plot == 'predom_supp':

            # set up supplemental predominance figure
            fig_predom_supp = plt.figure(figsize=(24, 13))
            ax = fig_predom_supp.add_subplot(1,1,1)

            map_predom(ax, predom,
                       xlim=None, ylim=None,
                       covars=covars,
                       which_covars='all',
                       cbar_axlab='predominant covar',
                       polys_color='none',
                       cbar_ticklab_fontsize=15,
                      )
            fig_predom_supp.subplots_adjust(left=0.02,
                                            right = 0.98,
                                            bottom=0.08,
                                            top=0.98,
                                            hspace=0.8,
                                            wspace=0.0
                                           )
            fig_predom_supp.savefig('FIG_SUPP_predom_all_covars.png', dpi=500)



    ##################
    # CREATE SUPP FIGS
    ##################

    #------------------
    # plot model errors

    if what_to_plot == 'error_supp':
        fig_err = plt.figure(figsize=(16,8))
        err_filename = 'err_map_yCOORDS_NIRv_100km.tif'
        rast = rxr.open_rasterio(os.path.join(data_dir,
                                              err_filename), masked=True)[0]
        # NOTE: multiply raster by -1 because I accidentally subtracted real value from
        #       prediction instead of vice versa
        ax = fig_err.add_subplot(111)
        divider = make_axes_locatable(ax)
        where = 'bottom'
        orientation = 'horizontal'
        size = '4%'
        cax = divider.append_axes(where, size=size, pad=0.2)

        vminmax = np.max(np.abs([np.nanpercentile(rast, 1),
                                 np.nanpercentile(rast, 99)]))
        rast.plot.imshow(ax=ax,
                         vmin=-vminmax,
                         vmax=vminmax,
                         cmap='cmo.balance_r',
                         add_colorbar=True,
                         cbar_ax=cax,
                         cbar_kwargs = {'orientation': orientation},
                         zorder=0,
                        )
        cax.tick_params(labelsize=cbarticklabel_fontsize)
        if orientation == 'vertical':
            axis = 'y'
        else:
            axis = 'x'
        getattr(cax, f'set_{axis}label')('standardized prediction error ($\Delta NIR_{V}/\Delta m$)', supp_axlabel_fontdict)
        phf.plot_juris_bounds(lev1_linewidth=0.1,
                              lev1_alpha=0.5,
                              lev1_zorder=1,
                              lev0_linewidth=0.25,
                              lev0_alpha=0.6,
                              lev0_zorder=2,
                              crs=plot_crs,
                              strip_axes=True,
                             )
        # format axes
        #ax.set_xlim(rast.rio.bounds()[0::2])
        #ax.set_ylim(rast.rio.bounds()[1::2])
        # NOTE: chopping off western edge because the equal earth projection
        #       makes NZ appear twice
        #ax.set_xlim(0.95 * ax.get_xlim()[0], ax.get_xlim()[1])

        # trim to the raster's bounds
        ax.set_xlim(rast.rio.bounds()[::2])
        ax.set_ylim(rast.rio.bounds()[1::2])

        # adjust subplots and save
        fig_err.subplots_adjust(bottom=0.08,
                                top=1,
                                left=0.02,
                                right=0.98,
                                wspace=0,
                                hspace=0,
                               )
        fig_err.savefig(f'FIG_SUPP_rf_error.png', dpi=600)


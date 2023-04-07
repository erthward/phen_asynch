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
from palettable.cartocolors.sequential import PinkYl_7
from palettable.scientific.sequential import Acton_20_r
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
#what_to_plot = 'main'
#what_to_plot = 'all_asynch'
what_to_plot = 'error_and_import'


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

# TODO: FIX EQUAL EARTH PROJECTION?
plot_crs = 8857

# load shapefiles
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                    'NewWorldFile_2020.shp')).to_crs(plot_crs)
# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                                  f)).to_crs(plot_crs))
subnational = pd.concat(subnational)


# asynch params for main results:
# main variable
var = 'NIRv'
# asynch neighborhood radius
neigh_rad = 100
# whether or not to use results of the model that included geocoods as covars
include_coords = 'y'

# SHAP-plotting params:
# whether or not to only use the top covars in the summary map
only_top_covars = True
# asynch percentile threshold below which to drop pixels from summary map
asynch_thresh = 75
# process the data raw, or read in processed files?
process_raw = False

# asynch colormap
asynch_cmap = 'cmo.dense_r'

# and predominance colormap
# NOTE: adapted from 'light' colormap at https://personal.sron.nl/~pault/#sec:qualitative
colors = [
          '#77AADD', # light blue
          '#FFAABB', # pink
          '#BBCC33', # pear
          '#EE8866', # orange
          '#EEDD88', # light yellow
          '#99DDFF', # light cyan
          '#44BB99', # mint
          #'#AAAA00', # olive
          #'#ADADAD', # light grey
         ]
if only_top_covars:
    colors = colors[:3]
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
top_covars = ['ppt.asy', 'tmp.min.asy', 'veg.ent']

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
    'no.predom': 'no predominant\ncovariate',
}

# stack all covars' SHAP maps
files = [f for f in os.listdir(data_dir) if re.search(
    f'SHAP_map_{include_coords}COORDS_[cdptv].*_{var}_{neigh_rad}km.tif', f)]
# reorder with top covars first, so that map colors of top covars are same
# color regardless of whether other colors are included in the map
reordered_files = []
for covar in covar_cbar_labels.keys():
    for f in files:
        if re.search(f'{covar}_{var}', f):
            reordered_files.append(f)
all_files = reordered_files
if only_top_covars:
    files = [f for f in all_files if re.search(f'(?<=COORDS_).*(?=_{var})', f).group() in top_covars]
else:
    files = all_files
shap_maps = []
covars = []
for f in files:
    covar = re.search(f'(?<=COORDS_).*(?=_{var})', f).group()
    covars.append(covar)
    shap_map = rxr.open_rasterio(os.path.join(data_dir, f), masked=True)[0]
    shap_map.name = covar
    shap_maps.append(shap_map)
xr_stack = xr.concat(shap_maps, dim=pd.Index(covars, name='covar'))
da = dask.array.stack(xr_stack)

if process_raw:
    # rechunk (to run slowly but successfully on laptop)
    # NOTE: axis sizes apparently needn't be exact multiples of their chunk sizes!
    da = da.rechunk(chunks=(da.shape[0], 1, 53, 500))

    # saturation: scaled standard deviation (as a metric of top covar predominance)
    # NOTE: 1 = max predominance; 0 = max evenness of importance
    std = deepcopy(shap_maps[0])
    std.values = minmax_scale(np.std(da, axis=0), 1, 99)
    assert np.nanmin(std) == 0 and np.nanmax(std) == 1
    # hue: predominant var
    predom = deepcopy(shap_maps[0])
    predom.values = np.argmax(np.abs(da), axis=0).compute()
    # NOTE: expressed in 360-degree values for hue (but cast as decimals,
    #       as required by the colorsys.hsv_to_rgb function I'm using), 
    #       yellow = 59/359, cyan = 179/359, and magenta = 299/359
    #predom = (59 + (120*predom))/359
    # mask to the std values
    predom = predom.where(pd.notnull(std))
    #assert np.all(np.unique(predom)[~np.isnan(np.unique(predom))] ==
    #              np.array([i/359 for i in [59,179,299]]))

else:
    if only_top_covars:
        file_suffix = 'top'
    else:
        file_suffix = 'all'
    std = rxr.open_rasterio(os.path.join(data_dir,
                    f'SHAP_std_{file_suffix}.tif'), masked=True)[0]
    predom = rxr.open_rasterio(os.path.join(data_dir,
                    f'SHAP_predom_{file_suffix}.tif'), masked=True)[0]

# reproject rasters
std = std.rio.write_crs(4326).rio.reproject(plot_crs, nodata=np.nan)
predom = predom.rio.write_crs(4326).rio.reproject(plot_crs, nodata=np.nan)

# asynchrony will be to mask out lower-asynchrony pixels
asynch = minmax_scale(rxr.open_rasterio(os.path.join(data_dir,
        f'{var}_STRICT_asynch_{neigh_rad}km.tif'), masked=True)[0], 0, 99)
# mask out pixels less than the percentile threshold
asynch = asynch.where(asynch >= np.nanpercentile(asynch, asynch_thresh))
# reproject and rescale
asynch = asynch.rio.reproject_match(std)
asynch = minmax_scale(asynch, 1, 99)
assert np.nanmin(asynch) == 0 and np.nanmax(asynch) == 1

# mask other two layers to the thresholded asynch map
std = std.where(pd.notnull(asynch))
predom = predom.where(pd.notnull(asynch))

# save rasters (if processed raw)
if process_raw:
    if only_top_covars:
        file_suffix = 'top'
    else:
        file_suffix = 'all'
    for lyr_name, lyr in {'std': std, 'predom': predom}.items():
        lyr.rio.to_raster(os.path.join(data_dir,
                                       f'SHAP_{lyr_name}_{file_suffix}.tif'))


# set up main figure
fig_main = plt.figure(figsize=(24, 26))
gs = fig_main.add_gridspec(ncols=270, nrows=270)


def plot_juris_bounds(ax, black_zorder=0, subnat_zorder=1, nat_zorder=2):
    countries.to_crs(plot_crs).plot(ax=ax,
                                color='black',
                                edgecolor='black',
                                linewidth=0,
                                alpha=0.6,
                                zorder=black_zorder,
                               )
    subnational.to_crs(plot_crs).plot(ax=ax,
                                  color='none',
                                  edgecolor='gray',
                                  linewidth=0.3,
                                  alpha=0.7,
                                  zorder=subnat_zorder,
                                 )
    countries.to_crs(plot_crs).plot(ax=ax,
                                color='none',
                                edgecolor='gray',
                                linewidth=0.5,
                                alpha=0.8,
                                zorder=nat_zorder,
                               )


def map_asynch(fig, cbar_axlab,
               gs=None, main_fig=True, var='NIRv',
               cbar_axlab_fontsize=cbar_axlab_fontsize,
               cbar_ticklab_fontsize=cbar_ticklab_fontsize):

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

        # add jurisdictional boundaries
        if main_fig:
            plot_juris_bounds(ax, 0, 2, 3)

        # read in the raster data and prepare it
        rast = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                              file), masked=True)[0]
        rast = rast.rio.write_crs(4326).rio.reproject(plot_crs)
                # NOTE: annoying AttributeError is because da.attrs['long_name']
        #       is retained as a tuple of names (rather than being subsetted
        #       by indexing) when I index a single layer out of an
        #       xarray.core.dataarray.DataArray;
        #       for now, a hacky fix is just assigning a string to that attr
        rast.attrs['long_name'] = ''
        rast.plot.imshow(ax=ax,
                         cmap=asynch_cmap,
                         vmin=np.nanpercentile(rast, 1),
                         vmax=np.nanpercentile(rast, 99),
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
        ax.set_xlim(rast.rio.bounds()[0::2])
        ax.set_ylim(rast.rio.bounds()[1::2])
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

        del rast


def map_rf_summ(ax, std, predom, xlim, ylim,
                cbar_axlab='covariate',
                cbar_axlab_fontsize=cbar_axlab_fontsize,
                cbar_ticklab_fontsize=cbar_ticklab_fontsize):
    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes('bottom', size='7%', pad=0.2)
    predom.plot.imshow(ax=ax,
                       cmap=predom_cmap,
                       add_colorbar=False,
                       #cbar_ax=cbar_ax,
                       #cbar_kwargs={'orientation': 'horizontal'},
                       zorder=1,
                      )
    ax.set_xticks(())
    ax.set_yticks(())

    for i in range(len(covars)):
        tmp_cmap = LinearSegmentedColormap.from_list('', ['#ffffff', colors[i]])
        std.where(predom==i).plot.imshow(ax=ax,
                                         cmap=tmp_cmap,
                                         add_colorbar=False,
                                         #vmin=0,
                                         #vmax=1,
                                         zorder=2+i
                                        )
    plot_juris_bounds(ax, 0, 5, 6)
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # make manual colorbar
    polys = []
    tick_locs = []
    tick_labs = []
    for i, covar in enumerate(covars+['no.predom']):
        polys.append(mplPolygon([(i, 0),
                              (i, 1),
                              (i+1, 1),
                              (i+1, 0),
                              (i, 0)]))
        tick_locs.append(i+0.5)
        tick_labs.append(covar_longnames[covar])
    patches = PatchCollection(polys, alpha=1, edgecolor='k')
    patches.set_color([mpl.colors.hex2color(c) for c in colors + ['#ffffff']])
    cbar_ax.add_collection(patches)
    cbar_ax.set_xlabel(cbar_axlab, fontdict={'fontsize': cbar_axlab_fontsize})
    cbar_ax.set_xticks(tick_locs)
    cbar_ax.set_xticklabels(tick_labs, fontdict={'fontsize':
                                                 cbar_ticklab_fontsize})
    cbar_ax.set_ylabel('')
    cbar_ax.set_yticks(())
    cbar_ax.set_xlim((0, len(covars)+1))
    cbar_ax.set_ylim((0,1))


cbar_axlab_dict = {'NIRv main': '$NIR_{V}\ asynchrony\ (\Delta NIR_{V}/\Delta m)$',
                   'NIRv': '$NIR_{V}\ asynch\ (\Delta NIR_{V}/\Delta m)$',
                   'SIF': '$SIF\  asynch\ (\Delta (mW m^{-2} sr^{-1} nm^{-1})/\Delta m)$',
                   'tmmn': '$tmp_{min}\ asynch\ (\Delta ^{\circ} C/\Delta m)$',
                   'tmmx': '$tmp_{max}\ asynch\ (\Delta ^{\circ} C/\Delta m)$',
                   'pr': '$ppt\ asynch\ (\Delta mm/\Delta m)$',
                   'def': '$cwd\ asynch\ (\Delta mm/\Delta m)$',
                   'cloud': '$cloud\ asynch\ (\Delta \%\ cover/\Delta m)$',

                  }

if __name__ == '__main__':

    #################
    # CREATE MAIN FIG
    #################

    plt.close('all')

    if what_to_plot == 'main':
        # plot the main asynch map
        map_asynch(fig_main, cbar_axlab_dict['NIRv main'], gs=gs, main_fig=True, var='NIRv')

        # plot the RF summary map
        ax = fig_main.add_subplot(gs[140:,:])
        map_rf_summ(ax, std, predom,
                    fig_main.axes[0].get_xlim(), fig_main.axes[0].get_ylim())

        # finish plot formatting
        fig_main.subplots_adjust(left=0.02,
                                 right = 0.98,
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

        fig_main.savefig('FIG_3_asynch_and_rf_results.png', dpi=500)

        del fig_main


    ##################
    # CREATE SUPP FIGS
    ##################

    #------------------------------
    # plot asynch maps for all vars

    first_supp_fig_num = 8
    vars = ['NIRv', 'SIF', 'tmmn', 'tmmx', 'pr', 'def', 'cloud']

    if what_to_plot == 'all_asynch':
        # make both vars' supp figs (each one stacking all 3 neighborhood radii)
        for n, var in enumerate(vars):
            print('\n\nNOW PRODUCING SUPPLEMENTAL FIG FOR %s..\n\n' % var)
            fig_supp = plt.figure(figsize=(9,12))
            map_asynch(fig_supp, cbar_axlab_dict[var],
                       gs=None, main_fig=False, var=var,
                       cbar_axlab_fontsize=13, cbar_ticklab_fontsize=10)
            fig_supp.subplots_adjust(bottom=0.02, top=0.95, left=0.02, right=0.88)
            fig_supp.savefig('FIG_S%i_%s_asynch_maps.png' % (first_supp_fig_num+n, var), dpi=600)
            del fig_supp



    #------------------
    # plot model errors

    if what_to_plot == 'error_and_import':
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
        subnational.to_crs(plot_crs).plot(ax=ax,
                                          color='none',
                                          edgecolor='black',
                                          linewidth=0.1,
                                          alpha=0.5,
                                          zorder=1,
                                         )
        countries.to_crs(plot_crs).plot(ax=ax,
                                        color='none',
                                        edgecolor='black',
                                        linewidth=0.25,
                                        alpha=0.6,
                                        zorder=2,
                                       )
        # format axes
        #ax.set_xlim(rast.rio.bounds()[0::2])
        #ax.set_ylim(rast.rio.bounds()[1::2])
        # NOTE: chopping off western edge because the equal earth projection
        #       makes NZ appear twice
        #ax.set_xlim(0.95 * ax.get_xlim()[0], ax.get_xlim()[1])
        ax.set_title('')
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_xticks(())
        ax.set_yticks(())

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
        sfig_num = first_supp_fig_num + len(vars) + 0
        fig_err.savefig(f'FIG_S{sfig_num}_rf_error.png', dpi=600)


        #-------------------------------
        # plot covar importance bar plots

        fig_import = plt.figure(figsize=(12,8))
        gs = fig_import.add_gridspec(ncols=4, nrows=3)
        for i, neigh_rad in enumerate([50, 100, 150]):
            for j, var in enumerate(['NIRv', 'SIF']):
                for k, include_coords in enumerate(['y', 'n']):
                    ax = fig_import.add_subplot(gs[i, j + (2*k)])
                    dfs = []
                    for import_metric in ['SHAP', 'permut']:
                        import_df = pd.read_csv(os.path.join(data_dir,
                            f'rf_{import_metric}_importance_{include_coords}COORDS_{var}_{neigh_rad}km.csv'))
                        import_df.columns = ['covariate', 'importance']
                        import_df.sort_values(by='importance', ascending=False, inplace=True)
                        import_df['metric'] = import_metric
                        # minmax normalize, to be able to display on common axis
                        minmax_scale = lambda vals: (vals-np.min(vals))/(np.max(vals)-np.min(vals))
                        import_df['importance'] = minmax_scale(import_df['importance'])
                        dfs.append(import_df)
                    df = pd.concat(dfs)
                    # make top rows the SHAP rows, sorted by descending importance
                    df = df.sort_values(by=['metric', 'importance'], ascending=[True, False])
                    g = sns.barplot(data=df,
                                    orient="h",
                                    x="importance",
                                    y="covariate",
                                    hue="metric",
                                    palette=['#6e6e6e', '#c4c4c4'],
                                    ci=None,
                                    alpha=1,
                                    ax=ax,
                                   )
                    if i==2 and j==1 and k==1:
                        pass
                    else:
                        ax.legend_.remove()
                    ax.set_ylabel('')
                    ax.set_xlabel('scaled importance' * int(i==2),
                                  fontdict=supp_axlabel_fontdict)
                    ax.tick_params(labelsize=cbarticklabel_fontsize)
                    if i == 0:
                        ax.set_title(f'{var}{"(coords)"*(include_coords=="y")}',
                                     fontdict=supp_axlabel_fontdict)
                    if j==0 and k==0:
                        ax.set_ylabel(f'{neigh_rad} km neigh. rad.',
                                      fontdict=supp_axlabel_fontdict)
                    #colors = [colorsys.hsv_to_rgb((59 + (129*i))/359, 1, 1) for i in range(3)]
                    # NOTE: reversing in order, to pop in correct oreder
                    #colors = colors[::-1]
                    for tick in ax.get_yticklabels():
                        tick.set_rotation(30)
                        tick.set_fontsize(cbarticklabel_fontsize)
                        if tick.get_text() in top_covars:
                            tick.set_weight('bold')

        # adjust subplots and save
        fig_import.subplots_adjust(bottom=0.10,
                                top=0.96,
                                left=0.14,
                                right=0.94,
                                wspace=0.6,
                                hspace=0.4,
                               )

        sfig_num = first_supp_fig_num + len(vars) + 1
        fig_import.savefig(f'FIG_S{sfig_num}_RF_var_import.png', dpi=400)



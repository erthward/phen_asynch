import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import geopandas as gpd
import rioxarray as rxr
from matplotlib.patches import Polygon as mplPolygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
from collections import Counter as C
import palettable
import cmocean
from palettable.cartocolors.sequential import PinkYl_7
from palettable.scientific.sequential import Acton_20_r
from shapely.geometry import Polygon as shapelyPolygon
from math import pi
from nlmpy import nlmpy
import statsmodels.api as sm
import sys
import re
import os

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf



# plot params
title_fontsize = 12
rowlab_fontsize = 18
axislab_fontsize = 11
ticklab_fontsize = 12
cbarlab_fontsize = 12
cbar_ticklab_fontsize = 10
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
arr_cmap = plt.cm.twilight_shifted
arr_mask_color = '#555555'
curve_cmap = plt.cm.bone
#asynch_cmap = 'viridis'
asynch_cmap = 'cmo.thermal'
#asynch_cmap = 'cmo.dense_r'
#asynch_cmap = mpl.cm.plasma.copy()
#asynch_cmap.set_bad('#717171')
rand_pix_to_track = [155, 67]
rand_pix_color = ['red', 'orange']
rand_pix_line_zorder = [40, 30]
rand_pix_linewidth = 1.5
#rand_pix_marker = '$\U0001F4CD$'
rand_pix_marker = '*'
rand_pix_markersize=80
scat_markersize=5
scat_markeralpha=0.75
#scat_label_text = {'high': '$steep\ =\ high\ asynch$',
#                  'low': '$flat\ =\ low\ asynch$',
#                 }
betas = [5, 1, 2, 1, 1]
noise_max = 0.25
rad = 10
min_x=0
max_x=1
min_y=0
max_y=1
dims=(21,21)
central_cell = [int((n-1)/2) for n in dims]
seed = 7031287
mpd_h = 1.4
orientation='landscape'
savefig=True


# load shapefiles
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                       'NewWorldFile_2020.shp'))
# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR,f)).to_crs(8857))
subnational = pd.concat(subnational)


# define functions
def get_seasonal_curve(betas, plot=True, color='gray',
                       ax=None, linestyle='-', linewidth=1, alpha=0.5,
                       min_seas_dist=None, max_seas_dist=None, zorder=None):
    if plot:
        if ax is None:
            no_ax = True
            fig, ax = plt.subplots()
        else:
            no_ax = False
    fitted = (betas[0] +
              betas[1]*np.sin(np.linspace(0,2*pi,365)) +
              betas[2]*np.cos(np.linspace(0,2*pi,365)) +
              betas[3]*np.sin(np.linspace(0,4*pi,365)) +
              betas[4]*np.cos(np.linspace(0,4*pi,365)))
    if plot:
        if isinstance(color, np.ndarray):
            seas_dist = pythag_dist(fitted, color)
            color_val = int(((seas_dist-min_seas_dist)/(
                            max_seas_dist-min_seas_dist)) * 255)
            #color=curve_cmap(color_val)
            color = 'black'
        ax.plot(np.linspace(1,365,365), fitted, linestyle,
                color=color, linewidth=linewidth, alpha=alpha, zorder=zorder)
        if no_ax:
            plt.show()

    return fitted


def make_betas_array(betas, asynch, dims=(21,21)):
    assert asynch in ['high', 'low']
    if asynch == 'high':
        arr = np.stack([nlmpy.mpd(*dims, h=mpd_h)-0.5*betas[i] for i in range(5)])
    else:
        stack = []
        for i in range(5):
            central_val = np.ones(dims)*betas[i]
            noise = np.random.uniform(0, noise_max, np.prod(dims)).reshape(dims)
            lyr = central_val + noise
            stack.append(lyr)
        arr = np.stack(stack)

    return(arr)


def pythag_dist(x, y):
    dist = np.sqrt(np.sum([(x[i]-y[i])**2 for i in range(len(x))]))
    return dist


def plot_all(betas, rad=rad, dims=(21,21), plot_it=True,
             min_seas_dist=0, max_seas_dist=1):
    # set the seed
    np.random.seed(seed)

    if plot_it:
        fig = plt.figure(dpi=dpi, figsize=(fig_width, fig_height))
        gs = fig.add_gridspec(nrows=4, ncols=4, width_ratios = [1.1, 0.2, 1, 1])
        top_axs = [fig.add_subplot(gs[0,i]) for i in [0,2,3]]
        bot_axs = [fig.add_subplot(gs[1,i]) for i in [0,2,3]]
    else:
        top_axs = [0,0]
        bot_axs = [0,0]
        max_seas_dists = []
        min_seas_dists = []
    for asynch, row_axs in zip(['low', 'high'], [top_axs, bot_axs]):
        if plot_it:
            ax1, ax2, ax3 = row_axs
        else:
            ax1 = ax2 = ax3 = None
        arr = make_betas_array(betas=betas, asynch=asynch, dims=dims)
           # array to gather values to plot in ax 1
        seas_dist_arr = np.ones(dims)*np.nan
        seas_peak_arr = np.ones(dims)*np.nan
        geo_dist_arr = np.ones(dims)*np.nan

        # lists to gather values to plot in ax 3
        xs = []
        ys = []
        cols = []

        cent_fitted = get_seasonal_curve(arr[:, central_cell[0], central_cell[1]],
                                             plot=False)
        pix_n = 0
        rand_pix_i = []
        rand_pix_j = []
        rand_pix_geo_dist = []
        rand_pix_seas_dist = []
        for i in range(arr.shape[1]):
            for j in range(arr.shape[2]):
                dist = pythag_dist([i,j], central_cell)
                geo_dist_arr[i,j] = dist
                linewidth = neighbor_linewidth
                alpha = neighbor_alpha
                #color_val = int(((rad-dist)/(rad-0)) * 255)
                #color=curve_cmap(color_val)
                fitted = get_seasonal_curve(arr[:, i, j], color=cent_fitted,
                                                ax=ax2, linewidth=linewidth,
                                                alpha=alpha,
                                            min_seas_dist=min_seas_dist,
                                            max_seas_dist=max_seas_dist,
                        # NOTE: only plot non-central cells and cells within rad
                        plot=([i, j] != central_cell and dist<rad and plot_it))

                # gather values for ax1 plot
                seas_dist = pythag_dist(fitted, cent_fitted)
                seas_dist_arr[i,j] = seas_dist
                seas_peak_arr[i,j] = np.argmax(fitted)

                # gather values for ax3 plot
                if dist <= rad:
                    xs.append(dist)
                    ys.append(seas_dist)
                    #cols.append(color)
                    cols.append(np.argmax(fitted))

                # save the random pixel to be tracked, if it's the right time
                if pix_n in rand_pix_to_track:
                    rand_pix_i.append(i)
                    rand_pix_j.append(j)
                    rand_pix_geo_dist.append(dist)
                    rand_pix_seas_dist.append(seas_dist)

                # increment the pixel counter
                pix_n += 1
        # now plot the central cell's curve
        cent_fitted = get_seasonal_curve(arr[:, central_cell[0], central_cell[1]],
                                         ax=ax2,
                                         color=central_curve_color,
                                         alpha=central_alpha,
                                         linewidth=central_linewidth,
                                         plot=plot_it)
        # and plot the random pixel curves
        for num, rand_pix_col in enumerate(rand_pix_color):
            fitted = get_seasonal_curve(arr[:, rand_pix_i[num], rand_pix_j[num]],
                                        color='k',
                                        ax=ax2, linewidth=1.5*rand_pix_linewidth,
                                        alpha=1, plot=plot_it,
                                        zorder=rand_pix_line_zorder[num])
            fitted = get_seasonal_curve(arr[:, rand_pix_i[num], rand_pix_j[num]],
                                        color=rand_pix_col,
                                        ax=ax2, linewidth=rand_pix_linewidth,
                                        alpha=1, plot=plot_it,
                                        zorder=rand_pix_line_zorder[num]+1)



        # plot image of seasonal distance and label with geo dist and radius)
        if plot_it:
            # ax1
            im = ax1.imshow(seas_peak_arr, cmap=arr_cmap,
                            vmin=0, vmax=364)
            cbar = plt.colorbar(im, ax=ax1, orientation='vertical')
            #cbar.set_label(label='time of year',
            cbar.set_label(label='',
                           size=cbarlab_fontsize,
                           rotation=-90,
                          )
            cbar.set_ticks(np.linspace(0,364,5))
            cbar.set_ticklabels(['Jan', 'Apr', 'Jul', 'Oct', 'Jan'])
            cbar.ax.tick_params(labelsize=cbar_ticklab_fontsize)
            cbar.ax.tick_params(labelsize=cbar_ticklab_fontsize, length=0)
            ax1.imshow(np.invert(geo_dist_arr>rad), cmap=plt.cm.gray,
                       vmin=0, vmax=1, alpha=rad_mask_alpha)
            #for i in range(dims[0]):
            #    for j in range(dims[1]):
            #        ax1.text(i, j, '%0.1f' % geo_dist_arr[i,j])
            rad_xs = [central_cell[1]+np.cos(ang)*rad for ang in np.linspace(0,
                                                                             2*pi, 720)]
            rad_ys = [central_cell[0]+np.sin(ang)*rad for ang in np.linspace(0,
                                                                             2*pi, 720)]
            ax1.plot(rad_xs, rad_ys, rad_linestyle,
                     color=rad_color, linewidth=rad_linewidth, alpha=rad_alpha)
            ax1.scatter(*central_cell, c='black', marker='*', s=rand_pix_markersize)
            for num, rand_pix_col in enumerate(rand_pix_color):
                ax1.scatter(rand_pix_i[num], rand_pix_j[num], c=rand_pix_col,
                            s=rand_pix_markersize, marker=rand_pix_marker,
                            edgecolor='k', linewidth=0.5)

            # ax3
            # NOTE: adjusting ys downward a bit, to make the low-high
            #       difference clearer
            reduction = 1
            ys = [y - reduction for y in ys]
            y = np.array(ys).T
            X = np.array(xs).T
            mod = sm.OLS(y, X, hasconst=False).fit()
            pred_xs = np.linspace(0.9*np.min(xs), 1.1*np.max(xs), 10)
            preds = mod.predict(pred_xs)
            ax3.scatter(xs, ys, c='k', s=scat_markersize,
                        alpha=scat_markeralpha)
            for num, rand_pix_col in enumerate(rand_pix_color):
                ax3.scatter(rand_pix_geo_dist[num], rand_pix_seas_dist[num]-reduction,
                            c=rand_pix_col, s=rand_pix_markersize,
                            marker=rand_pix_marker, edgecolor='k', linewidth=0.5)
            # cover the (0,0) point, to avoid confusion
            ax3.scatter(0, 0, c='white', s=2*scat_markersize, edgecolor='white')
            if asynch == 'low':
                ax3.set_ylim((min_seas_dist, max_seas_dist))
            else:
                ax3.set_ylim((min_seas_dist, max_seas_dist))
            ax3.plot(pred_xs, preds, '-k')
            # make sure origin is bottom-left corner
            ax3.set_xlim(0, ax3.get_xlim()[1])

            # dress it up
            ax1.set_xticks(())
            ax1.set_yticks(())
            ax1.set_xticklabels(())
            ax1.set_yticklabels(())
            ax1.tick_params(axis=u'both', which=u'both', length=0)
            if asynch == 'high':
                ax1.set_xlabel('date of peak phenology',
                               fontdict={'fontsize': axislab_fontsize})
            ax1.set_ylabel('%s\nasynchrony' % asynch,
                           fontdict={'fontsize': rowlab_fontsize})

            ax2.set_xticks(())
            ax2.set_xticklabels(())
            ax2.set_yticks(())
            ax2.set_yticklabels(())
            if asynch == 'high':
                ax2.set_xlabel('time of year',
                           fontdict={'fontsize': axislab_fontsize})
            ax2.set_ylabel('phenological signal',
                           fontdict={'fontsize': axislab_fontsize})

            ax3.set_xticks(())
            ax3.set_xticklabels(())
            ax3.set_yticks(())
            ax3.set_yticklabels(())
            if asynch == 'high':
                ax3.set_xlabel('geographic distance',
                          fontdict={'fontsize': axislab_fontsize})
            ax3.set_ylabel('phenological distance',
                           fontdict={'fontsize': axislab_fontsize})

        # if not plotting, just store the min and max values
        else:
            max_seas_dists.append(np.max(ys))
            min_seas_dists.append(np.min(ys))

    if not plot_it:
        return (np.min(min_seas_dists), np.max(max_seas_dists))
    else:
        # adjust suplot spacing
        plt.subplots_adjust(left=subplots_adj_left,
                        bottom=subplots_adj_bottom,
                        right=subplots_adj_right,
                        top=subplots_adj_top,
                        wspace=subplots_adj_wspace,
                        hspace=subplots_adj_hspace)
        return(fig, gs, mod)


def map_asynch(fig, cbar_axlab,
               gs=None, main_fig=True, var='NIRv',
               cbar_axlab_fontsize=18, cbar_ticklab_fontsize=10):

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

        # either grab the lower half of the main fig
        if main_fig:
            #ax = fig.add_subplot(gs[3:, :])
            ax = fig.add_subplot(2,1,2)
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

        # read in the raster data and prepare it
        rast = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                              file), masked=True)[0]
        rast = rast.rio.write_crs(4326).rio.reproject(8857)
                # NOTE: annoying AttributeError is because da.attrs['long_name']
        #       is retained as a tuple of names (rather than being subsetted
        #       by indexing) when I index a single layer out of an
        #       xarray.core.dataarray.DataArray;
        #       for now, a hacky fix is just assigning a string to that attr
        rast.attrs['long_name'] = ''
        rast.plot.imshow(ax=ax,
                         zorder=0,
                         cmap=asynch_cmap,
                         vmin=np.nanpercentile(rast, 1),
                         vmax=np.nanpercentile(rast, 99),
                         add_colorbar=True,
                         cbar_ax=cax,
                         cbar_kwargs = {'orientation': orientation},
                        )
        cax.tick_params(labelsize=cbar_ticklab_fontsize)
        if main_fig:
            cax.set_xlabel(cbar_axlab, fontdict={'fontsize': cbar_axlab_fontsize})
        else:
            cax.set_ylabel(cbar_axlab, fontdict={'fontsize': cbar_axlab_fontsize})
        subnational.to_crs(8857).plot(ax=ax,
                                      color='none',
                                      edgecolor='black',
                                      zorder=1,
                                      alpha=0.6,
                                     )
        countries.to_crs(8857).plot(ax=ax,
                                    color='none',
                                    edgecolor='black',
                                    linewidth=1,
                                    alpha=0.8,
                                    zorder=2,
                                   )
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
                         fontdict={'fontsize': 42})

        del rast

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
    plt.close('all')
    # make the conceptual figure
    print('\n\nNOW PRODUCING FIG 3...\n\n')
    min_seas_dist, max_seas_dist = plot_all(betas, rad=rad, dims=dims,
                                            plot_it=False)
    fig, gs, mod = plot_all(betas, rad=rad, dims=dims, min_seas_dist=min_seas_dist,
                   max_seas_dist=max_seas_dist, plot_it=True)
    # add the asynch map below
    map_asynch(fig, cbar_axlab_dict['NIRv main'], gs=gs, main_fig=True, var='NIRv')

    # add labels for parts A. and B.
    fig.axes[0].text(-5.8, -5, 'A.', size=24, weight='bold')
    fig.axes[-2].text(1.11*fig.axes[-2].get_xlim()[0],
                      1.065*fig.axes[-2].get_ylim()[1],
                      'B.', size=24, weight='bold')

    # adjust subplots and save
    fig.subplots_adjust(bottom=0.05, top=0.92, left=0.04, right=0.98)
    fig.savefig('FIG_3_asynch_concept_and_map.png', dpi=600)

    # make both vars' supp figs (each one stacking all 3 neighborhood radii)
    for n, var in enumerate(['NIRv', 'SIF', 'tmmn', 'tmmx', 'pr', 'def', 'cloud']):
        print('\n\nNOW PRODUCING SUPPLEMENTAL FIG FOR %s..\n\n' % var)
        fig_supp = plt.figure(figsize=(19,24))
        map_asynch(fig_supp, cbar_axlab_dict[var],
                   gs=None, main_fig=False, var=var,
                   cbar_axlab_fontsize=26, cbar_ticklab_fontsize=20)
        fig_supp.subplots_adjust(bottom=0.02, top=0.95, left=0.0, right=0.88)
        fig_supp.savefig('FIG_S%i_%s_asynch_maps.png' % (6+n, var), dpi=600)


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
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf


# plot params
title_fontsize = 11
rowlab_fontsize = 14
axislab_fontsize = 12
ticklab_fontsize = 12
cbarlab_fontsize = 12
cbar_ticklab_fontsize = 10
annot_fontsize = 11
scat_label_fontsize = 8
scat_label_fontcolor = 'black'
scat_label_linewidth = 3
fig_width = 10.5
fig_height = 4
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
rand_pix_to_track = [155, 89]
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
        gs = fig.add_gridspec(nrows=20, ncols=80)
        top_axs = [fig.add_subplot(gs[top_pt:bot_pt, l_pt:r_pt]) for
                   top_pt, bot_pt, l_pt, r_pt in [[0, 7, 8, 21],
                                                  [0, 9, 25, 51],
                                                  [0, 9, 55, 95]]]
        bot_axs = [fig.add_subplot(gs[top_pt:bot_pt, l_pt:r_pt]) for
                   top_pt, bot_pt, l_pt, r_pt in [[13, 20, 8, 21],
                                                  [11, 20, 25, 51],
                                                  [11, 20, 55, 95]]]

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
                ax1.set_xlabel('')
            ax1.set_ylabel((f'{" "*0*(asynch=="high")}{asynch}'
                            f'\nasynchrony{" "*0*(asynch=="low")}'),
                           labelpad=70,
                           rotation=0,
                           fontdict={'fontsize': rowlab_fontsize,
                                     'weight': 'bold'},
                          )

            ax2.set_xticks(())
            ax2.set_xticklabels(())
            ax2.set_yticks(())
            ax2.set_yticklabels(())
            if asynch == 'high':
                ax2.set_xlabel('time of year',
                           fontdict={'fontsize': axislab_fontsize})
            if asynch=='high':
                ax2.set_ylabel(f'{" "*35}phenological signal',
                               fontdict={'fontsize': axislab_fontsize})

            ax3.set_xticks(())
            ax3.set_xticklabels(())
            ax3.set_yticks(())
            ax3.set_yticklabels(())
            if asynch == 'high':
                ax3.set_xlabel('geographic distance',
                          fontdict={'fontsize': axislab_fontsize})
                ax3.set_ylabel(f'{" "*35}phenological distance',
                           fontdict={'fontsize': axislab_fontsize})

        # if not plotting, just store the min and max values
        else:
            max_seas_dists.append(np.max(ys))
            min_seas_dists.append(np.min(ys))

    # add circular colorbar
    if plot_it and asynch=='high':
        cbar_ax = fig.add_subplot(gs[7:13, 12:17], projection='polar')
        azimuths = np.arange(90, 451, 1)
        zeniths = np.arange(40, 70, 1)
        values = azimuths * np.ones((30, 361))
        cbar_ax.pcolormesh(azimuths*np.pi/180.0, zeniths, values,
                           cmap=mpl.cm.twilight_shifted_r)
        cbar_ax.set_xticks(np.array([90, 180, 270, 360])/180*np.pi,
                           ['Jan', 'Oct', 'Jul', 'Apr'],
                           size=6,
                          )
        cbar_ax.tick_params(pad=0.1)
        cbar_ax.set_yticks(())
        ax1.text(-5,
                 7,
                 'date of peak phenology ',
                 size=axislab_fontsize,
                 clip_on=False,
                 rotation=90,
                )
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


if __name__ == '__main__':
    plt.close('all')
    # make the conceptual figure
    print('\n\nNOW PRODUCING CONCEPTUAL ASYNCH FIG...\n\n')
    min_seas_dist, max_seas_dist = plot_all(betas, rad=rad, dims=dims,
                                            plot_it=False)
    fig, gs, mod = plot_all(betas, rad=rad, dims=dims, min_seas_dist=min_seas_dist,
                   max_seas_dist=max_seas_dist, plot_it=True)

    # adjust subplots and save
    fig.subplots_adjust(bottom=0.06, top=0.92, left=0.04, right=0.98)
    fig.savefig(os.path.join(phf.FIGS_DIR, 'FIG_SUPP_asynch_concept.png'), dpi=600)


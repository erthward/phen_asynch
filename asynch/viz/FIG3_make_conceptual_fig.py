import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as mplPolygon
from collections import Counter as C
from palettable.cartocolors.sequential import PinkYl_7
from palettable.scientific.sequential import Acton_20_r
from shapely.geometry import Polygon as shapelyPolygon
from math import pi
from nlmpy import nlmpy
import statsmodels.api as sm
import sys
import re
import os


# plot params
title_fontsize = 12
rowlab_fontsize = 18
axislab_fontsize = 11
ticklab_fontsize = 12
cbarlab_fontsize = 12
cbarticklab_fontsize = 10
annot_fontsize = 14
scat_label_fontsize = 10
scat_label_fontcolor = 'black'
scat_label_linewidth = 3
fig_width = 10.5
fig_height = 5.6
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
neighbor_linewidth = 1
neighbor_alpha = 0.5
rad_linestyle = ':'
rad_color = 'white'
rad_linewidth = 2
rad_alpha = 0.75
rad_mask_alpha = 0.175
arr_cmap = plt.cm.twilight_shifted
arr_mask_color = '#555555'
curve_cmap = plt.cm.viridis
asynch_scat_cmap = curve_cmap
rand_pix_to_track = 208
rand_pix_color = 'red'
rand_pix_linewidth = 1.5
#rand_pix_marker = '$\U0001F4CD$'
rand_pix_marker = '*'
rand_pix_markersize=80
scat_markersize=5
scat_markeralpha=0.75
scat_label_text = {'high': '$steep\ =\ high\ asynch$',
                   'low': '$flat\ =\ low\ asynch$',
                  }
betas = [5, 1, 2, 1, 1]
noise_max = 0.25
rad = 10
min_x=0
max_x=1
min_y=0
max_y=1
dims=(21,21)
central_cell = [int((n-1)/2) for n in dims]
seed = 703121887
mpd_h = 1.4
orientation='landscape'
savefig=True

def get_seasonal_curve(betas, plot=True, color='gray',
                       ax=None, linestyle='-', linewidth=1, alpha=0.5,
                       min_seas_dist=None, max_seas_dist=None):
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
            color=curve_cmap(color_val)
        ax.plot(np.linspace(1,365,365), fitted, linestyle,
                color=color, linewidth=linewidth, alpha=alpha)
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
        gs = fig.add_gridspec(nrows=2, ncols=4, width_ratios = [1.1, 0.2, 1, 1])
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
                if pix_n == rand_pix_to_track:
                    rand_pix_i = i
                    rand_pix_j = j
                    rand_pix_geo_dist = dist
                    rand_pix_seas_dist = seas_dist

                # increment the pixel counter
                pix_n += 1
        # now plot the central cell's curve
        cent_fitted = get_seasonal_curve(arr[:, central_cell[0], central_cell[1]],
                                         ax=ax2,
                                         color=central_curve_color,
                                         alpha=central_alpha,
                                         linewidth=central_linewidth,
                                         plot=plot_it)
        # and plot the random pixel's curve
        fitted = get_seasonal_curve(arr[:, rand_pix_i, rand_pix_j],
                                    color='k',
                                    ax=ax2, linewidth=1.5*rand_pix_linewidth,
                                    alpha=1, plot=plot_it)
        fitted = get_seasonal_curve(arr[:, rand_pix_i, rand_pix_j],
                                    color=rand_pix_color,
                                    ax=ax2, linewidth=rand_pix_linewidth,
                                    alpha=1, plot=plot_it)



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
            cbar.ax.tick_params(labelsize=cbarticklab_fontsize)
            cbar.ax.tick_params(labelsize=cbarticklab_fontsize, length=0)
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
            ax1.scatter(*central_cell, c='black', s=rand_pix_markersize)
            ax1.scatter(rand_pix_i, rand_pix_j, c=rand_pix_color,
                        s=rand_pix_markersize, marker=rand_pix_marker,
                        edgecolor='k', linewidth=0.5)

            # ax3
            y = np.array(ys).T
            X = np.stack((np.ones((len(xs))), np.array(xs))).T
            #X = np.array(xs).T
            mod = sm.OLS(y, X, hasconst=True).fit()
            #mod = sm.OLS(y, X, hasconst=False).fit()
            pred_xs = np.stack((np.ones(10), np.linspace(0.9*np.min(xs), 1.1*np.max(xs), 10))).T
            #pred_xs = np.linspace(0.9*np.min(xs), 1.1*np.max(xs), 10)
            preds = mod.predict(pred_xs)
            #ax3.scatter(xs, ys, c=cols, s=scat_markersize,
            #            cmap=asynch_scat_cmap)
            ax3.scatter(xs, ys, c='k', s=scat_markersize,
                        alpha=scat_markeralpha)
            ax3.scatter(rand_pix_geo_dist, rand_pix_seas_dist,
                        c=rand_pix_color, s=rand_pix_markersize,
                        marker=rand_pix_marker, edgecolor='k', linewidth=0.5)
            # cover the (0,0) point, to avoid confusion
            ax3.scatter(0, 0, c='white', s=2*scat_markersize, edgecolor='white')
            if asynch == 'low':
                ax3.set_ylim((min_seas_dist, max_seas_dist))
            else:
                ax3.set_ylim((min_seas_dist, max_seas_dist))
            ax3.plot(pred_xs[:,1], preds, '-', color=scat_label_fontcolor)
            #ax3.plot(pred_xs, preds, '-k')
            #ax3.plot([pred_xs[2,1], pred_xs[2,1], pred_xs[5,1]],
            #ax3.plot([pred_xs[2], pred_xs[2], pred_xs[5]],
            #         [preds[2], preds[5], preds[5]], ':',
            #         color=scat_label_fontcolor,
            #         linewidth=scat_label_linewidth)
            ax3.text(0.05, 0.85, scat_label_text[asynch],
                     transform=ax3.transAxes,
                     fontdict={'fontsize': scat_label_fontsize,
                               'color': scat_label_fontcolor})


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
        return(fig, mod)


if __name__ == '__main__':
    plt.close('all')
    min_seas_dist, max_seas_dist = plot_all(betas, rad=rad, dims=dims,
                                            plot_it=False)
    fig, mod = plot_all(betas, rad=rad, dims=dims, min_seas_dist=min_seas_dist,
                   max_seas_dist=max_seas_dist, plot_it=True)
    if savefig:
        fig.savefig('seasonal_asynch_conceptual_figure.png',
                    dpi=dpi, orientation=orientation)
    else:
        fig.show()


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
rowlab_fontsize = 20
axislab_fontsize = 11
ticklab_fontsize = 12
annot_fontsize = 14
scat_label_fontsize = 6
scat_label_fontcolor = 'red'
fig_width = 10.5
fig_height = 5.6
dpi = 400
n_ticklabels = 5
subplots_adj_left=0.05
subplots_adj_bottom=0.05
subplots_adj_right=0.96
subplots_adj_top=0.95
subplots_adj_wspace=0.14
subplots_adj_hspace=0.30
central_curve_color = '#220000'
central_linewidth = 3
central_alpha = 1
neighbor_linewidth = 1
neighbor_alpha = 0.5
rad_linestyle = ':'
rad_color = 'white'
rad_linewidth = 2
rad_alpha = 0.75
rad_mask_alpha = 0.175
arr_cmap = plt.cm.viridis
arr_mask_color = '#555555'
curve_cmap = plt.cm.viridis_r
scat_label_text = {'high': 'steep slope:\nhigh asynch',
                   'low': 'shallow slope:\nlow asynch',
                  }
betas = [1, 2, 1, 1]
noise_max = 0.3
rad = 10
min_x=0
max_x=1
min_y=0
max_y=1
dims=(21,21)
central_cell = [int((n-1)/2) for n in dims]
seed = 9
mpd_h = 1.2
orientation='landscape'
savefig=True

def get_seasonal_curve(betas, plot=True, color='gray',
                       ax=None, linestyle='-', linewidth=1, alpha=0.5):
    if plot:
        if ax is None:
            no_ax = True
            fig, ax = plt.subplots()
        else:
            no_ax = False
    fitted = (betas[0]*np.sin(np.linspace(0,2*pi,365)) +
              betas[1]*np.cos(np.linspace(0,2*pi,365)) +
              betas[2]*np.sin(np.linspace(0,4*pi,365)) +
              betas[3]*np.cos(np.linspace(0,4*pi,365)))
    if plot:
        ax.plot(np.linspace(1,365,365), fitted, linestyle,
                color=color, linewidth=linewidth, alpha=alpha)
        if no_ax:
            plt.show()

    return fitted


def make_betas_array(betas, asynch, dims=(21,21)):
    assert asynch in ['high', 'low']
    if asynch == 'high':
        arr = np.stack([nlmpy.mpd(*dims, h=mpd_h)-0.5*betas[i] for i in range(4)])
    else:
        stack = []
        for i in range(4):
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
        top_axs = [fig.add_subplot(2,3,i+1) for i in range(3)]
        bot_axs = [fig.add_subplot(2,3,i+4) for i in range(3)]
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
        for i in range(arr.shape[1]):
            for j in range(arr.shape[2]):
                dist = pythag_dist([i,j], central_cell)
                geo_dist_arr[i,j] = dist
                linewidth = neighbor_linewidth
                alpha = neighbor_alpha
                color_val = int(((rad-dist)/(rad-0)) * 255)
                color=curve_cmap(color_val)
                fitted = get_seasonal_curve(arr[:, i, j], color=color,
                                                ax=ax2, linewidth=linewidth,
                                                alpha=alpha,
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
                    cols.append(color)
        # now plot the central cell's curve
        cent_fitted = get_seasonal_curve(arr[:, central_cell[0], central_cell[1]],
                                         ax=ax2,
                                         color=central_curve_color,
                                         alpha=central_alpha,
                                         linewidth=central_linewidth,
                                         plot=plot_it)

        # plot image of seasonal distance and label with geo dist and radius)
        if plot_it:
            ax1.imshow(seas_dist_arr, cmap=arr_cmap,
                       vmin=min_seas_dist, vmax=max_seas_dist)
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


            y = np.array(ys).T
            X = np.array(xs).T
            mod = sm.OLS(y, X).fit()
            pred_xs = np.linspace(0.9*np.min(xs), 1.1*np.max(xs), 10)
            preds = mod.predict(pred_xs)
            ax3.scatter(xs, ys, c=cols)
            ax3.set_ylim((min_seas_dist, max_seas_dist))
            ax3.plot(pred_xs, preds, '-k')
            ax3.plot([pred_xs[4], pred_xs[7], pred_xs[7]],
                     [preds[4], preds[4], preds[7]], ':r')
            ax3.text(pred_xs[7], preds[2], scat_label_text[asynch],
                     fontdict={'fontsize': scat_label_fontsize,
                               'color': scat_label_fontcolor})


            # dress it up
            ax1.set_xticks(())
            ax1.set_yticks(())
            ax1.set_xticklabels(())
            ax1.set_yticklabels(())
            ax1.tick_params(axis=u'both', which=u'both', length=0)
            if asynch == 'high':
                ax1.set_xlabel('map of seasonal distance',
                               fontdict={'fontsize': axislab_fontsize})

            # dress it up
            ax1.set_xticks(())
            ax1.set_yticks(())
            ax1.set_xticklabels(())
            ax1.set_yticklabels(())
            ax1.tick_params(axis=u'both', which=u'both', length=0)
            if asynch == 'high':
                ax1.set_xlabel('map of seasonal peak',
                               fontdict={'fontsize': axislab_fontsize})
            ax1.set_ylabel('%s asynch' % asynch,
                           fontdict={'fontsize': rowlab_fontsize})

            ax2.set_xticks(())
            ax2.set_xticklabels(())
            ax2.set_yticks(())
            ax2.set_yticklabels(())
            if asynch == 'high':
                ax2.set_xlabel('day of year',
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


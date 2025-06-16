#!/bin/python

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os, sys

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                                        'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf


def make_corr_plot(corr_dict,
                   fig,
                   ax,
                   linewidth=5,
                   cmap='cividis',
                   add_cbar=False,
                   cbar_ax=None,
                   neigh_text_fontsize=7,
                   r2_text_fontsize=6,
                   cbar_axlab_fontsize=9,
                   cbar_ticklab_fontsize=5,
                  ):
    '''
    plot equilteral triangle of correlation values betweenasynchrony
    neighborhood distances (100 km always being at top and 50 km and 150 km
    being the left and right bases)
    '''
    # validate corr_dict
    assert len(corr_dict) == 3
    corr_dict_sort = {tuple(np.sort(k).ravel()): v for k, v in corr_dict.items()}
    comps = (50, 100), (100, 150), (50, 150)
    for comp in comps:
        assert comp in corr_dict_sort.keys()

    # plot
    neigh_positions = {100: (0.5, 1.1),
                       50: (-0.1, 0),
                       150: (1.1, 0),
                      }
    line_positions = {(50, 100):  [(0.1, 0.1), (0.4, 0.9)],
                      (100, 150): [(0.9, 0.1), (0.6, 0.9)],
                      (50, 150):  [(0.1, 0.0), (0.9, 0.0)],
                     }
    r2_positions = {(50, 100): (0.0, 0.5),
                    (100, 150): (1.0, 0.5),
                    (50, 150): (0.50, -0.25),
                   }
    for neigh, pos in neigh_positions.items():
        ax.text(*pos,
                f"{neigh}\nkm",
                size=neigh_text_fontsize,
                horizontalalignment='center',
                verticalalignment='center',
               )
    line_segs = []
    corrs = []
    colors = []
    for comp in comps:
        line_seg = line_positions[comp]
        line_segs.append(line_seg)
        corrs.append(corr_dict_sort[comp])
        colors.append(mpl.colormaps.get(cmap)(corr_dict_sort[comp]))
    lc = LineCollection(line_segs,
                        linewidths=[linewidth]*len(line_segs),
                       )

    lc.set_colors(colors)
    lc.set_array(np.asarray(corrs))
    lc.set_cmap(cmap)
    ax.add_collection(lc)
    if add_cbar:
        cmap = getattr(mpl.cm, cmap)
        norm = mpl.colors.Normalize(vmin=0, vmax=1)
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        assert cbar_ax is not None
        cbar = plt.colorbar(sm,
                            cax=cbar_ax,
                            orientation='horizontal',
                           )
        cbar.ax.set_xlabel('$R^2$', size=cbar_axlab_fontsize, loc='left')
        cbar.ax.set_xticks([0, 0.25, 0.5, 0.75, 1],
                           ['0', '0.25', '0.50', '0.75', '1'],
                          )
        cbar.ax.tick_params(labelsize=cbar_ticklab_fontsize)
    for comp, r2 in corr_dict_sort.items():
        ax.text(*r2_positions[comp],
                "%0.2f" % r2,
                size=r2_text_fontsize,
                horizontalalignment='center',
                verticalalignment='center',
               )
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_xlim(-0.15, 1.15)
    ax.set_ylim(-0.25, 1.55)
    ax.set_axis_off()


# load inter-neighborhood correlations for all variables
corr_df = pd.read_csv(os.path.join(phf.TABS_DIR,
                                   "TAB_asynch_neigh_rad_comp_r2s.csv",
                                  ))

# create figure and fill it with correlation triangle plots
fig = plt.figure(figsize=(1, 10))
gs = fig.add_gridspec(ncols=1,
                      nrows=len(np.unique(corr_df['var']))+2,
                      height_ratios=[0.05, 0.1] + [1]*len(np.unique(corr_df['var'])),
                     )
cbar_ax=fig.add_subplot(gs[0,:])
for i, var in enumerate(corr_df['var'].unique()):
    sub_df = corr_df[corr_df['var'] == var]
    str_corr_dict = dict([*zip(sub_df['neigh_comp'].values,
                               sub_df['r2'].values)])
    corr_dict = {tuple([int(i) for i in k.split('|')]): v for k, v in
                 str_corr_dict.items()}
    ax = fig.add_subplot(gs[i+2,:])
    make_corr_plot(corr_dict=corr_dict,
                   fig=fig,
                   ax=ax,
                   add_cbar=i==(len(np.unique(corr_df['var']))-1),
                   cbar_ax=cbar_ax,
                  )
fig.subplots_adjust(bottom=0.01,
                    top=0.99,
                    left=0.08,
                    right=0.92,
                    hspace=0,
                    wspace=0,
                   )
fig.savefig(os.path.join(phf.FIGS_DIR,
                         'FIG_SUPP_asynch_map_neigh_corr_triangles.png'),
            dpi=500,
           )


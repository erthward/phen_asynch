import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from copy import deepcopy
import rioxarray as rxr
import seaborn as sns
import cmocean
import dask
import os, sys, re

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf



"""
TODO:
    - clip all maps to the same bounding box used (the one in other figs, if
            equal area)
    - save, then add, other importance metric
    - add letters to label parts
    - other TODOs
    - if keeping interp_map, could like colors in the map to font colors or
    something used to label rows at the right
    - in interp map, use asynch-clustering from final analysis to mask map instead of
      pixel-wise asynch values?

                            IMPORT          VAR         SHAP

              PRED MAP

                            IMPORT          VAR         SHAP



                                            VAR          SHAP



                                            VAR          SHAP
                    SHAP INTERP


                                            VAR          SHAP





"""


# set plotting params
# TODO: FIX EQUAL EARTH PROJECTION?
#plot_crs = 8857
plot_crs = 4326

title_fontdict = {'fontsize': 30}
maps_axlabel_fontdict = {'fontsize': 23}
other_axlabel_fontdict = {'fontsize': 20}
cbarlabel_fontdict = {'fontsize': 14}
ticklabel_fontsize=12
cbarticklabel_fontsize=8

# set paths
data_dir = phf.EXTERNAL_RF_DATA_DIR

# load shapefiles
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                       'NewWorldFile_2020.shp'))
# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                                  f)).to_crs(plot_crs))
subnational = pd.concat(subnational)

var = 'NIRv'
neigh_rad = 100

top_covars = {'tmp.min.nsd': ('CHELSA_bio6_1981-2010_V.2.1_5km_'
                              '10CELLRAD_NEIGHSD.tif'),
              'ppt.sea.nsd': ('CHELSA_bio15_1981-2010_V.2.1_5km_'
                              '10CELLRAD_NEIGHSD.tif'),
              'ppt.asy': f'pr_asynch_{neigh_rad}km.tif',
              'tmp.min.asy': f'tmmn_asynch_{neigh_rad}km.tif',
             }
#top_covar_ax_labels = {
#    'tmp.min.nsd': '$s_{ngh} T_{min,ann}$',
#    'ppt.sea.nsd': '$s_{ngh} P_{ann}$',
#    'ppt.asy':     '$asynch_{P}$',
#    'tmp.min.asy': '$asynch_{T_{min,mon}}$'
#                       }
top_covar_cbar_labels = {
    'tmp.min.nsd': '$(^{\circ}C)$',
    'ppt.sea.nsd': '$(mm)$',
    'ppt.asy': '$\Delta dist_{seas_{P}}/\Delta  dist_{geo}$',
    'tmp.min.asy': '$\Delta dist_{seas_{T_{min}}}/\Delta  dist_{geo}$',
                        }

# set up figure
fig = plt.figure(figsize=(26.25, 15.25))
gs = fig.add_gridspec(ncols=260, nrows=150)


#------------------------------------------
# 1. plot original vars and their SHAP maps
# NOTE: currently hard-coded for 4 top covars
row_idxs = np.array([0, 35, 70, 105])
row_start_idxs = row_idxs + 10
row_end_idxs = row_idxs + 10 + 30
col_start_idxs = [140, 205]
col_end_idxs = [190, 255]

# loop over all SHAP tifs, to determine the min and max SHAP vals for plotting
vmins = []
vmaxs = []
for i, top_covar_item in enumerate(top_covars.items()):
    top_covar, covar_filename = top_covar_item
    shap_filename = f'SHAP_map_{top_covar}_{var}_{neigh_rad}km.tif'
    rast = rxr.open_rasterio(os.path.join(data_dir,
                                          shap_filename), masked=True)[0]
    vmin = np.nanpercentile(rast, 1)
    vmax = np.nanpercentile(rast, 99)
    vmins.append(vmin)
    vmaxs.append(vmax)
vminmax = np.max(np.abs(vmins+vmaxs))
shap_vmin = -vminmax
shap_vmax = vminmax


def change_cbar_ticklabels_to_scinot(cax, orientation):
    if orientation == 'horizontal':
        axis = 'x'
    else:
        axis = 'y'
    ticks = getattr(cax, f'get_{axis}ticks')()
    new_ticklabs = ['%0.1e' % t for t in ticks]
    getattr(cax, f'set_{axis}ticks')(ticks, new_ticklabs)


for i, top_covar_item in enumerate(top_covars.items()):
    top_covar, covar_filename = top_covar_item
    shap_filename = f'SHAP_map_{top_covar}_{var}_{neigh_rad}km.tif'
    ax_covar = fig.add_subplot(gs[row_start_idxs[i]:row_end_idxs[i],
                                  col_start_idxs[0]:col_end_idxs[0]])
    ax_shap = fig.add_subplot(gs[row_start_idxs[i]:row_end_idxs[i],
                                 col_start_idxs[1]:col_end_idxs[1]])
    cmaps = ['cmo.matter_r', 'cmo.curl_r']
    j = 0
    for ax, filename in zip([ax_covar, ax_shap],
                            [covar_filename, shap_filename]):
        divider = make_axes_locatable(ax)
        where = 'bottom'
        orientation = 'horizontal'
        size = '4%'
        cax = divider.append_axes(where, size=size, pad=0.2)
        rast = rxr.open_rasterio(os.path.join(data_dir,
                                              filename), masked=True)[0]
        rast = rast.rio.write_crs(4326).rio.reproject(plot_crs)
        # NOTE: hack to get around the stupid xarray AttributeError
        rast.attrs['long_name'] = ''
        if j == 1:
            vmin = shap_vmin
            vmax = shap_vmax
        else:
            vmin = np.nanpercentile(rast, 1)
            vmax = np.nanpercentile(rast, 99)
        rast.plot.imshow(ax=ax,
                         zorder=0,
                         cmap=cmaps[j],
                         vmin=vmin,
                         vmax=vmax,
                         add_colorbar=True,
                         cbar_ax=cax,
                         cbar_kwargs = {'orientation': orientation},
                        )
        cax.tick_params(labelsize=cbarticklabel_fontsize)
        if j == 0:
            cbarlabel = top_covar_cbar_labels[top_covar]
        else:
            cbarlabel = ''
        if orientation == 'vertical':
            cax.set_ylabel(cbarlabel, cbarlabel_fontdict)
        else:
            cax.set_xlabel(cbarlabel, cbarlabel_fontdict)
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
        ax.set_xlim(rast.rio.bounds()[0::2])
        ax.set_ylim(rast.rio.bounds()[1::2])
        # NOTE: chopping off western edge because the equal earth projection
        #       makes NZ appear twice
        ax.set_xlim(0.95 * ax.get_xlim()[0], ax.get_xlim()[1])
        if i == 0:
            if j == 0:
                ax.set_title('covariate', fontdict=maps_axlabel_fontdict)
            else:
                ax.set_title('SHAP values', fontdict=maps_axlabel_fontdict)
        else:
            ax.set_title('')
        if j == 0:
            #ax.set_ylabel(top_covar_ax_labels[top_covar],
            ax.set_ylabel(top_covar, fontdict=maps_axlabel_fontdict)
        else:
            ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_xticks(())
        ax.set_yticks(())
        del rast
        j += 1




#----------------------------------------------
# 2. TODO: PLOT PREDS - ACTUAL VALS IN TOP LEFT

ax = fig.add_subplot(gs[10:55, 0:90])
ax.text(1,1, "PREDS-VALS HERE!", size=20)
ax.set_xlim(0, 5)
ax.set_ylim(0,5)




#-----------------------------------
# 3. plot covar importance bar plots

dfs = []
# TODO: ADD IN OTHER IMPORTANCE METRIC
#for import_metric in ['SHAP', 'info']:
for import_metric in ['SHAP', 'SHAP']:
    import_df = pd.read_csv(os.path.join(data_dir,
                    f'rf_{import_metric}_importance_{var}_{neigh_rad}km.csv'))
    import_df.columns = ['covariate', 'importance']
    import_df.sort_values(by='importance', ascending=False, inplace=True)
    import_df['metric'] = import_metric
    dfs.append(import_df)
    # TODO: NEED TO STANDARDIZE EACH METRIC SO THAT THEY CAN SHARE THE X?
df = pd.concat(dfs)
# make top rows the SHAP rows, sorted by descending importance
df = df.sort_values(by=['metric', 'importance'], ascending=[True, False])

ax = fig.add_subplot(gs[10:60, 105:130])
g = sns.barplot(data=df,
                orient="h",
                x="importance",
                y="covariate",
                hue="metric",
                palette=['red', 'blue'],
                alpha=1,
                ax=ax,
               )
g.legend(loc='upper center',
         bbox_to_anchor=(0.25, 1.5),
         ncol=2,
         fontsize=ticklabel_fontsize,
         title='feature importance',
         title_fontsize=maps_axlabel_fontdict['fontsize'],
        )
ax.set_ylabel('')
ax.set_xlabel('importance', fontdict=cbarlabel_fontdict)
ax.tick_params(labelsize=ticklabel_fontsize)
for tick in ax.get_yticklabels():
    tick.set_rotation(30)
    if tick.get_text() in top_covars:
        tick.set_weight('bold')


#    ax.barh(df['covariate'],
#            df['importance'],
#            color=
#                orient='h',
#                palette='OrRd',
#                kwargs={'width': 2},
#                ax=ax,
#               )
#    ax.legend([], [], frameon=False)
#    ax.set_xlabel('RF covariate', fontdict=axlabel_fontdict)
#    ax.set_ylabel('importance', fontdict=axlabel_fontdict)
#    ax.set_title('')
#    ax.tick_params(labelsize=ticklabel_fontsize)
#    #for tick in ax.get_yticklabels():
#    #     tick.set_rotation(55)

    # TODO:
    # add labels for parts A. and B.
    #fig.axes[0].text(-5.8, -5, 'A.', size=24, weight='bold')
    #fig.axes[-2].text(1.11*fig.axes[-2].get_xlim()[0],
    #                 1.065*fig.axes[-2].get_ylim()[1],
    #                 'B.', size=24, weight='bold')


#--------------------------------
# 4. plot SHAP interpretation map

ax = fig.add_subplot(gs[70:, :130])
divider = make_axes_locatable(ax)
where = 'bottom'
orientation = 'horizontal'
size = '4%'
cax = divider.append_axes(where, size=size, pad=0.2)

# include lat and lon (i.e., y and x) maps?
include_latlon = False

# only use top-importance covars?
use_only_top_import_covars = False

# get SHAP summary map filename
filename = 'SHAP_interp_map_%s_%ikm%s%s.tif' % (var, neigh_rad,
                                                '_LATLON' * include_latlon,
                                    '_TOPIMPORT' * use_only_top_import_covars)
rast = rxr.open_rasterio(os.path.join(data_dir, filename), masked=True)[0]

# use a special color for the background?
bg_color = None
#bg_color = '#ffffff'

# set colormap for the covars
# NOTE: using 12-color palette at http://tsitsul.in/blog/coloropt/
color_hex = ['#ebac23',
             '#b80058',
             '#008cf9',
             '#006e00',
             '#00bbad',
             '#d163e6',
             '#b24502',
             '#ff9287',
             '#5954d6',
             '#00c6f8',
             '#878500',
             '#00a76c',
             '#bdbdbd',
            ]
cmap = ListedColormap(color_hex[:len(covar_names)])
# make missing values black (to make colors more discernable)
if bg_color is not None:
    cmap.set_bad(bg_color)

# plot it
countries.to_crs(rast.rio.crs).plot(color='black',
                                      ax=ax,
                                      zorder=0)
im = rast.plot.imshow(cmap=cmap,
                        ax=ax,
                        add_colorbar=True,
                        cbar_ax=cax,
                        zorder=1,
                       )
cbar = im.colorbar
ymin, ymax = cbar.ax.get_ylim()
tick_locs = np.linspace(ymax/len(covar_names)/2,
                        ymax - (ymax/len(covar_names)/2),
                        len(covar_names))
cbar.ax.set_yticks(tick_locs, covar_names)
# format
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('predomnant covars, by |SHAP value|', fontdict=title_fontdict)

# adjust subplots and save
fig.subplots_adjust(bottom=0.03,
                    top=0.92,
                    left=0.08,
                    right=0.9,
                    wspace=0,
                    hspace=0,
                   )

#fig.savefig('FIG_4_RF_summary.png', dpi=700)

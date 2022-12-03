import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as path_effects
from copy import deepcopy
import rioxarray as rxr
import seaborn as sns
import colorsys
import cmocean
import dask
import os, sys, re

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf



# indicate the main variable and neighborhood radius,
# and whether or not to use the model that included the geo-coord polynomial
var = 'NIRv'
neigh_rad = 100
include_coords = 'y'

# set plotting params
# TODO: FIX EQUAL EARTH PROJECTION?
#plot_crs = 8857
plot_crs = 4326

maps_axlabel_fontdict = {'fontsize': 30}
cbarlabel_fontdict = {'fontsize': 22}
cbarticklabel_fontsize = 18
partlabel_fontsize=40

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

top_covars = {'ppt.asy': f'pr_asynch_{neigh_rad}km.tif',
              'tmp.min.asy': f'tmmn_asynch_{neigh_rad}km.tif',
              'veg.ent': 'MODIS_IGBP_veg_entropy.tif',
             }
top_covar_cbar_labels = {
    'ppt.asy': '$\Delta mm/\Delta m$',
    'tmp.min.asy': '$\Delta ^{\circ} c/\Delta m$',
    'veg.ent': '$entropy$',
                        }

# set up figure
fig = plt.figure(figsize=(24, 18))
gs = fig.add_gridspec(ncols=260, nrows=170)



#---------------------
# 1. plot model errors
fig_err = plt.figure(figsize=(16,8))
err_filename = f'err_map_{include_coords}COORDS_{var}_{neigh_rad}km.tif'
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
getattr(cax, f'set_{axis}label')('standardized prediction error ($\Delta NIR_{V}/\Delta m$)', maps_axlabel_fontdict)
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

# save this Axes object's bounds, to use for other axes later on
common_xlim = ax.get_xlim()
common_ylim = ax.get_ylim()

# add part label
#ax.text(1.11 * ax.get_xlim()[0],
#        1.17 * ax.get_ylim()[1],
#        'A.', size=partlabel_fontsize, weight='bold')

# adjust subplots and save
fig_err.subplots_adjust(bottom=0.08,
                        top=1,
                        left=0.02,
                        right=0.98,
                        wspace=0,
                        hspace=0,
                       )
fig_err.savefig('FIG_S12_RF_err_map.png', dpi=700)



#-----------------------------------
# 2. plot covar importance bar plots

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

ax = fig.add_subplot(gs[10:70, 5:40])
g = sns.barplot(data=df,
                orient="h",
                x="importance",
                y="covariate",
                hue="metric",
                palette=['#6e6e6e', '#c4c4c4'],
                alpha=1,
                ax=ax,
               )
g.legend(loc='upper center',
         bbox_to_anchor=(0.5, 1.14),
         ncol=2,
         fontsize=cbarticklabel_fontsize,
         #title='feature importance',
         #title_fontsize=maps_axlabel_fontdict['fontsize'],
        )
ax.set_ylabel('')
ax.set_xlabel('scaled importance', fontdict=maps_axlabel_fontdict)
ax.tick_params(labelsize=cbarticklabel_fontsize)
#colors = [colorsys.hsv_to_rgb((59 + (129*i))/359, 1, 1) for i in range(3)]
# NOTE: reversing in order, to pop in correct oreder
#colors = colors[::-1]
for i, tick in enumerate(ax.get_yticklabels()):
    tick.set_rotation(30)
    tick.set_fontsize(cbarticklabel_fontsize)
    if tick.get_text() in top_covars:
        tick.set_weight('bold')
        # color to match the colors in part D's HSV map
        #color = colors.pop()
        #tick.set_path_effects([path_effects.Stroke(linewidth=3,
        #                                           foreground=color),
        #                       path_effects.Normal()])

# add part label
ax.text(-0.6, -1.06, 'A.', size=partlabel_fontsize, weight='bold')


#--------------------------------
# 3. plot SHAP interpretation map

ax = fig.add_subplot(gs[80:, :])
divider = make_axes_locatable(ax)
where = 'bottom'
orientation = 'horizontal'
size = '3%'
cax = divider.append_axes(where, size=size, pad=0.35)

# save to GeoTIFF
var_order = '_'.join([*top_covars])
hsv = rxr.open_rasterio(os.path.join(phf.EXTERNAL_RF_DATA_DIR,
    f'SHAP_hsv_map_{var_order}_{include_coords}COORDS_{var}_{neigh_rad}km.tif'),
                        masked=True)

# plot it
# NOTE: create black background
cmap = mpl.colors.LinearSegmentedColormap.from_list('black', ['#000000']*2)
cmap.set_bad('#000000')
black_bg = deepcopy(hsv[0])
black_bg.plot.imshow(ax=ax,
                     cmap=cmap,
                     vmin=0,
                     vmax=0,
                     add_colorbar=False,
                     zorder=0)
hsv.plot.imshow(ax=ax, zorder=1)
subnational.to_crs(hsv.rio.crs).plot(ax=ax,
                                      color='none',
                                      edgecolor='gray',
                                      linewidth=0.4,
                                      alpha=0.5,
                                      zorder=2,
                                     )
countries.to_crs(hsv.rio.crs).plot(ax=ax,
                                   color='none',
                                   edgecolor='gray',
                                   linewidth=0.6,
                                   alpha=0.8,
                                   zorder=3,
                                  )

if orientation == 'horizontal':
    axis = 'x'
else:
    axis = 'y'
axmin, axmax = getattr(cax, f'get_{axis}lim')()
# add 'no predominant covariate' to the covar names that will be used to label
# the colorbar
covar_names = [*top_covars]
covar_names.append('no predom.')
tick_locs = np.linspace(axmax/len(covar_names)/2,
                        axmax - (axmax/len(covar_names)/2),
                        len(covar_names))
if orientation == 'vertical':
    getattr(cax, f'set_{axis}ticks')(tick_locs[::-1], covar_names)
elif orientation == 'horizontal':
    getattr(cax, f'set_{axis}ticks')(tick_locs[::-1], covar_names[::-1])

# format
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('')
ax.set_xlim(hsv.rio.bounds()[::2])
ax.set_ylim(hsv.rio.bounds()[1::2])
cax.tick_params(length=0)
for i, tick in enumerate(getattr(cax, f'get_{axis}ticklabels')()):
    if orientation == 'vertical':
        tick.set_rotation(-45)
    tick.set_fontsize(cbarlabel_fontdict['fontsize'])
if orientation == 'horizontal':
    cax.set_xlabel('predominant covariate',
                   labelpad=1.1,
                   fontdict={**maps_axlabel_fontdict, 'rotation':0},
                  )
elif orientation == 'vertical':
    cax.set_ylabel('predominant\ncovariate',
                   labelpad=0.5,
                   fontdict={**maps_axlabel_fontdict, 'rotation':90},
                  )
if orientation == 'vertical':
    cax.set_xticks(())
elif orientation == 'horizontal':
    cax.set_yticks(())

# add custom colorbar patches
cbar_xmin, cbar_xmax = cax.get_xlim()
cbar_ymin, cbar_ymax = cax.get_ylim()
if orientation == 'horizontal':
    breaks = np.linspace(cbar_xmin, cbar_xmax, 5)
else:
    breaks = np.linspace(cbar_ymin, cbar_ymax, 5)
# create an RGB array containing yellow, cyan, magenta, and white
colors = np.array([colorsys.hsv_to_rgb((59 + (129*i))/359,
                                       1-(i==3), 1) for i in range(4)])
assert np.all(colors.shape == np.array((4,3)))
patches = []
for i in [*range(4)][::-1]:
    if orientation == 'horizontal':
        poly = Polygon([[breaks[::-1][i], cbar_ymin],
                        [breaks[::-1][i], cbar_ymax],
                        [breaks[::-1][i+1], cbar_ymax],
                        [breaks[::-1][i+1], cbar_ymin],
                        [breaks[::-1][i], cbar_ymin]])
    else:
        poly = Polygon([[cbar_xmin, breaks[i]],
                        [cbar_xmax, breaks[i]],
                        [cbar_xmax, breaks[i+1]],
                        [cbar_xmin, breaks[i+1]],
                        [cbar_xmin, breaks[i]]])
    patches.append(poly)
pc = PatchCollection(patches, alpha=1, edgecolor='k')
pc.set_color(colors)
cax.add_collection(pc)

# move colorbar axis and tick labels to the right
cax.yaxis.set_label_position("right")
cax.yaxis.tick_right()

# add part label
ax.text(1.15 * ax.get_xlim()[0],
        0.95 * ax.get_ylim()[1],
        'D.', size=partlabel_fontsize, weight='bold')



#------------------------------------------
# 4. plot original vars and their SHAP maps
# NOTE: CURRENTLY HARD-CODED FOR 3 VARS
col_idxs = np.array([55, 125, 195])
col_start_idxs = col_idxs + 5
col_end_idxs = col_idxs + 5 + 60
row_start_idxs = [5, 45]
row_end_idxs = [35, 75]

# loop over all SHAP tifs, to determine the min and max SHAP vals for plotting
vmins = []
vmaxs = []
for j, top_covar_item in enumerate(top_covars.items()):
    top_covar, covar_filename = top_covar_item
    shap_filename = f'SHAP_map_{include_coords}COORDS_{top_covar}_{var}_{neigh_rad}km.tif'
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

#colors = [colorsys.hsv_to_rgb((59 + (129*i))/359, 1, 1) for i in range(3)]
# NOTE: reversing in order, to pop in correct oreder
#colors = colors[::-1]
for j, top_covar_item in enumerate(top_covars.items()):
    top_covar, covar_filename = top_covar_item
    shap_filename = f'SHAP_map_{include_coords}COORDS_{top_covar}_{var}_{neigh_rad}km.tif'
    ax_covar = fig.add_subplot(gs[row_start_idxs[0]:row_end_idxs[0],
                                  col_start_idxs[j]:col_end_idxs[j]])
    ax_shap = fig.add_subplot(gs[row_start_idxs[1]:row_end_idxs[1],
                                 col_start_idxs[j]:col_end_idxs[j]])
    # NOTE: cmo.thermal matches the asynchrony map in fig 3
    cmaps = ['cmo.thermal', 'cmo.curl_r']
    i = 0
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
        # NOTE: get rid of egregiously large values about longitudinal bounds
        #       of Equal Earth projection
        #rast = rast.where(np.abs(rast<=1e37))
        # NOTE: hack to get around the stupid xarray AttributeError
        rast.attrs['long_name'] = ''
        if i == 1:
            vmin = shap_vmin
            vmax = shap_vmax
        else:
            vmin = np.nanpercentile(rast, 1)
            vmax = np.nanpercentile(rast, 99)
        rast.plot.imshow(ax=ax,
                         zorder=0,
                         cmap=cmaps[i],
                         vmin=vmin,
                         vmax=vmax,
                         add_colorbar=True,
                         cbar_ax=cax,
                         cbar_kwargs = {'orientation': orientation},
                        )
        cax.tick_params(labelsize=cbarticklabel_fontsize)
        if i == 0:
            cbarlabel = top_covar_cbar_labels[top_covar]
        else:
            cbarlabel = 'SHAP val'
        if orientation == 'vertical':
            axis = 'y'
        else:
            axis = 'x'
        getattr(cax, f'set_{axis}label')(cbarlabel, cbarlabel_fontdict)
        # convert cbar-axes ticklabels to scientific notation, if necessary
        #cax_ticklabs = cax.get_xticklabels()
        #if True in [len(tl)>=5 for tl in cax_ticklabs]:
        #    new_ticklabs = ['%0.2e' % (float(n)) for n in cax_ticklabs]
        #    cax.set_xticks(cax.get_xticks(), new_ticklabs)
        #cax.ticklabel_format(style='sci')
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
        ax.set_xlim(common_xlim)
        ax.set_ylim(common_ylim)
        if j == 0:
            if i == 0:
                ax.set_ylabel('covariate', fontdict=maps_axlabel_fontdict)
            else:
                ax.set_ylabel('influence', fontdict=maps_axlabel_fontdict)
        else:
            ax.set_ylabel('')
        if i == 0:
            #ax.set_ylabel(top_covar_ax_labels[top_covar],
            lab = ax.set_title(top_covar, fontdict=maps_axlabel_fontdict)
            #color = colors.pop()
            #lab.set_path_effects([path_effects.Stroke(linewidth=3,
            #                                       foreground=color),
            #                   path_effects.Normal()])
        else:
            ax.set_title('')
        ax.set_xlabel('')
        ax.set_xticks(())
        ax.set_yticks(())
        del rast

        # add part label
        if i == 0 and j == 0:

            ax.text(1.35 * ax.get_xlim()[0],
                    1.4 * ax.get_ylim()[1],
                    'B.', size=partlabel_fontsize, weight='bold')

        elif i == 1 and j == 0:

            ax.text(1.35 * ax.get_xlim()[0],
                    1.4 * ax.get_ylim()[1],
                    'C.', size=partlabel_fontsize, weight='bold')

        # set number of x-axis ticks on cbar axes
        if i == 0:
            nticks = 4
        else:
            nticks = 5
        cax.xaxis.set_major_locator(plt.MaxNLocator(nticks))

        i += 1


# adjust subplots and save
fig.subplots_adjust(bottom=0.02,
                    top=0.99,
                    left=0.06,
                    right=0.92,
                    wspace=0.02,
                    hspace=0.03,
                   )
print('saving...')
fig.savefig('FIG_4_RF_summary.png', dpi=600)


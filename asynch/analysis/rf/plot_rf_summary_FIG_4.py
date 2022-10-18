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

title_fontdict = {'fontsize': 30}
maps_axlabel_fontdict = {'fontsize': 20}
other_axlabel_fontdict = {'fontsize': 20}
cbarlabel_fontdict = {'fontsize': 14}
ticklabel_fontsize = 12
interp_map_ticklabel_fontsize = 14
cbarticklabel_fontsize = 8
partlabel_fontsize=30

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
    'ppt.asy': '$\Delta dist_{seas_{P}}/\Delta  dist_{geo}$',
    'tmp.min.asy': '$\Delta dist_{seas_{T_{min}}}/\Delta  dist_{geo}$',
    'veg.ent': '$entropy$',
                        }

# set up figure
fig = plt.figure(figsize=(26.25, 10.5))
gs = fig.add_gridspec(ncols=260, nrows=150)



#---------------------
# 1. plot model errors

err_filename = f'err_map_{include_coords}COORDS_{var}_{neigh_rad}km.tif'
rast = rxr.open_rasterio(os.path.join(data_dir,
                                      err_filename), masked=True)[0]
# NOTE: multiply raster by -1 because I accidentally subtracted real value from
#       prediction instead of vice versa
ax = fig.add_subplot(gs[10:55, 0:90])
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
getattr(cax, f'set_{axis}label')('prediction error', maps_axlabel_fontdict)
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
ax.text(1.11 * ax.get_xlim()[0],
        1.17 * ax.get_ylim()[1],
        'A.', size=partlabel_fontsize, weight='bold')



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

ax = fig.add_subplot(gs[10:50, 105:130])
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
         fontsize=11,
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
    tick.set_fontsize(10)
    if tick.get_text() in top_covars:
        tick.set_weight('bold')
        # color to match the colors in part D's HSV map
        #color = colors.pop()
        #tick.set_path_effects([path_effects.Stroke(linewidth=3,
        #                                           foreground=color),
        #                       path_effects.Normal()])

# add part label
ax.text(-0.35, -1.1, 'B.', size=partlabel_fontsize, weight='bold')


#--------------------------------
# 3. plot SHAP interpretation map

ax = fig.add_subplot(gs[60:140, :130])
divider = make_axes_locatable(ax)
where = 'bottom'
orientation = 'horizontal'
size = '7%'
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

if orientation == 'vertical':
    axis = 'y'
else:
    axis = 'x'
axmin, axmax = getattr(cax, f'get_{axis}lim')()
# add 'no predominant covariate' to the covar names that will be used to label
# the colorbar
covar_names = [*top_covars]
covar_names.append('no predom.')
tick_locs = np.linspace(axmax/len(covar_names)/2,
                        axmax - (axmax/len(covar_names)/2),
                        len(covar_names))
getattr(cax, f'set_{axis}ticks')(tick_locs, covar_names)

# format
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('')
ax.set_xlim(hsv.rio.bounds()[::2])
ax.set_ylim(hsv.rio.bounds()[1::2])
cax.tick_params(length=0, labelsize=interp_map_ticklabel_fontsize)
for i, tick in enumerate(cax.get_yticklabels()):
    tick.set_rotation(45)
    if i == 3:
        tick.set_style('italic')
cax.set_xlabel('predominant covariate', fontdict=maps_axlabel_fontdict)
cax.set_yticks(())

# add custom colorbar patches
cbar_xmin, cbar_xmax = cax.get_xlim()
cbar_ymin, cbar_ymax = cax.get_ylim()
breaks = np.linspace(cbar_xmin, cbar_xmax, 5)
# create an RGB array containing yellow, cyan, magenta, and white
colors = np.array([colorsys.hsv_to_rgb((59 + (129*i))/359,
                                       1-(i==3), 1) for i in range(4)])
assert np.all(colors.shape == np.array((4,3)))
patches = []
for i in range(4):
    poly = Polygon([[breaks[i], cbar_ymin],
                    [breaks[i], cbar_ymax],
                    [breaks[i+1], cbar_ymax],
                    [breaks[i+1], cbar_ymin],
                    [breaks[i], cbar_ymin]])
    patches.append(poly)
pc = PatchCollection(patches, alpha=1, edgecolor='k')
pc.set_color(colors)
cax.add_collection(pc)

# add part label
ax.text(1.06 * ax.get_xlim()[0],
        1.05 * ax.get_ylim()[1],
        'D.', size=partlabel_fontsize, weight='bold')



#------------------------------------------
# 4. plot original vars and their SHAP maps
# NOTE: CURRENTLY HARD-CODED FOR 3 VARS
row_idxs = np.array([0, 45, 90])
row_start_idxs = row_idxs + 5
row_end_idxs = row_idxs + 5 + 40
col_start_idxs = [150, 210]
col_end_idxs = [200, 260]

# loop over all SHAP tifs, to determine the min and max SHAP vals for plotting
vmins = []
vmaxs = []
for i, top_covar_item in enumerate(top_covars.items()):
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
for i, top_covar_item in enumerate(top_covars.items()):
    top_covar, covar_filename = top_covar_item
    shap_filename = f'SHAP_map_{include_coords}COORDS_{top_covar}_{var}_{neigh_rad}km.tif'
    ax_covar = fig.add_subplot(gs[row_start_idxs[i]:row_end_idxs[i],
                                  col_start_idxs[0]:col_end_idxs[0]])
    ax_shap = fig.add_subplot(gs[row_start_idxs[i]:row_end_idxs[i],
                                 col_start_idxs[1]:col_end_idxs[1]])
    # NOTE: cmo.thermal matches the asynchrony map in fig 3
    cmaps = ['cmo.thermal', 'cmo.curl_r']
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
            axis = 'y'
        else:
            axis = 'x'
        getattr(cax, f'set_{axis}label')(cbarlabel, cbarlabel_fontdict)
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
        if i == 0:
            if j == 0:
                ax.set_title('covariate', fontdict=maps_axlabel_fontdict)
            else:
                ax.set_title('SHAP values', fontdict=maps_axlabel_fontdict)
        else:
            ax.set_title('')
        if j == 0:
            #ax.set_ylabel(top_covar_ax_labels[top_covar],
            lab = ax.set_ylabel(top_covar, fontdict=maps_axlabel_fontdict)
            #color = colors.pop()
            #lab.set_path_effects([path_effects.Stroke(linewidth=3,
            #                                       foreground=color),
            #                   path_effects.Normal()])
        else:
            ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_xticks(())
        ax.set_yticks(())
        del rast

        # add part label
        if i == 0 and j == 0:

            ax.text(1.5 * ax.get_xlim()[0],
                    1.27 * ax.get_ylim()[1],
                    'C.', size=partlabel_fontsize, weight='bold')

        j += 1


# adjust subplots and save
fig.subplots_adjust(bottom=0,
                    top=1,
                    left=0.03,
                    right=0.98,
                    wspace=0,
                    hspace=0,
                   )

fig.savefig('FIG_4_RF_summary.png', dpi=700)


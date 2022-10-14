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
    - figure out equal-area projection
    - other TODOs
    - if keeping interp_map, could link colors in the map to font colors or
    something used to label rows at the right?
"""

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



#---------------------
# 1. plot model errors

err_filename = f'err_map_{include_coords}COORDS_{var}_{neigh_rad}km.tif'
rast = rxr.open_rasterio(os.path.join(data_dir,
                                      err_filename), masked=True)[0]
# NOTE: multiply raster by -1 because I accidentally subtracted real value from
#       prediction instead of vice versa
ax = fig.add_subplot(gs[10:55, 0:90])
divider = make_axes_locatable(ax)
where = 'top'
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
cax.xaxis.tick_top()
cax.xaxis.set_label_position('top')
if orientation == 'vertical':
    axis = 'y'
else:
    axis = 'x'
getattr(cax, f'set_{axis}label')('prediction error', cbarlabel_fontdict)
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

# save this Axes object's bounds, to use for other axes later on
common_xlim = ax.get_xlim()
common_ylim = ax.get_ylim()

# add part label
ax.text(1.11 * ax.get_xlim()[0],
        1.37 * ax.get_ylim()[1],
        'A.', size=partlabel_fontsize, weight='bold')



#-----------------------------------
# 2. plot covar importance bar plots

dfs = []
# TODO: ADD IN OTHER IMPORTANCE METRIC
#for import_metric in ['SHAP', 'info']:
for import_metric in ['SHAP', 'SHAP']:
    import_df = pd.read_csv(os.path.join(data_dir,
        f'rf_{import_metric}_importance_{include_coords}COORDS_{var}_{neigh_rad}km.csv'))
    import_df.columns = ['covariate', 'importance']
    import_df.sort_values(by='importance', ascending=False, inplace=True)
    import_df['metric'] = import_metric
    dfs.append(import_df)
    # TODO: NEED TO STANDARDIZE EACH METRIC SO THAT THEY CAN SHARE THE X?
df = pd.concat(dfs)
# make top rows the SHAP rows, sorted by descending importance
df = df.sort_values(by=['metric', 'importance'], ascending=[True, False])

ax = fig.add_subplot(gs[10:50, 105:130])
g = sns.barplot(data=df,
                orient="h",
                x="importance",
                y="covariate",
                hue="metric",
                palette=['#9fd6c5', '#9fafd6'],
                alpha=1,
                ax=ax,
               )
g.legend(loc='upper center',
         bbox_to_anchor=(0.5, 1.25),
         ncol=2,
         fontsize=cbarticklabel_fontsize,
         #title='feature importance',
         #title_fontsize=maps_axlabel_fontdict['fontsize'],
        )
ax.set_ylabel('')
ax.set_xlabel('importance', fontdict=cbarlabel_fontdict)
ax.tick_params(labelsize=cbarticklabel_fontsize)
for tick in ax.get_yticklabels():
    tick.set_rotation(30)
    tick.set_fontsize(10)
    if tick.get_text() in top_covars:
        tick.set_weight('bold')

# add part label
ax.text(-9000, -0.6, 'B.', size=partlabel_fontsize, weight='bold')


#--------------------------------
# 3. plot SHAP interpretation map

ax = fig.add_subplot(gs[60:140, :130])
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


# get the covar names (identically to how they were determined in the file that
# produced the map)
shap_files = [f for f in os.listdir(data_dir) if (f.startswith('SHAP_map')
                                                  and f.endswith('.tif'))]
shap_files = [f for f in shap_files if (var in f and
                                        str(neigh_rad)+'km' in f and
                                        f'{include_coords}COORDS' in f)]
if not include_latlon:
    shap_files = [f for f in shap_files if not
                      re.search('(?<=^SHAP_map_)[xy].?1?_%s' % var, f)]
if use_only_top_import_covars:
    shap_files = [f for f in shap_files if ('ppt' in f or
                                            'tmp.min' in f or
                        # NOTE: need x and y in case include_latlon == True
                            re.search('(?<=^SHAP_map_)[xy].?1?_%s' % var, f))]
# read all covars
covar_names = [re.search('(?<=^SHAP_map_).*(?=_%s)' % var,
                         f).group() for f in shap_files]


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
                      cbar_kwargs = {'orientation': orientation},
                      zorder=1,
                     )
cbar = im.colorbar
if orientation == 'vertical':
    axis = 'y'
else:
    axis = 'x'
axmin, axmax = getattr(cbar.ax, f'get_{axis}lim')()
tick_locs = np.linspace(axmax/len(covar_names)/2,
                        axmax - (axmax/len(covar_names)/2),
                        len(covar_names))
getattr(cbar.ax, f'set_{axis}ticks')(tick_locs, covar_names)
# format
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('')
cax.tick_params(length=0, labelsize=interp_map_ticklabel_fontsize)
for tick in cax.get_yticklabels():
    tick.set_rotation(45)

# add part label
ax.text(1.06 * ax.get_xlim()[0],
        1.05 * ax.get_ylim()[1],
        'D.', size=partlabel_fontsize, weight='bold')



#------------------------------------------
# 4. plot original vars and their SHAP maps
# NOTE: currently hard-coded for 4 top covars
row_idxs = np.array([0, 35, 70, 105])
row_start_idxs = row_idxs + 10
row_end_idxs = row_idxs + 10 + 30
col_start_idxs = [150, 210]
col_end_idxs = [200, 260]

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
            ax.set_ylabel(top_covar, fontdict=maps_axlabel_fontdict)
        else:
            ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_xticks(())
        ax.set_yticks(())
        del rast

        # add part label
        if i == 0 and j == 0:

            ax.text(1.5 * ax.get_xlim()[0],
                    1.35 * ax.get_ylim()[1],
                    'C.', size=partlabel_fontsize, weight='bold')

        j += 1


# adjust subplots and save
fig.subplots_adjust(bottom=0.01,
                    top=1,
                    left=0.03,
                    right=0.98,
                    wspace=0,
                    hspace=0,
                   )

fig.savefig('FIG_4_RF_summary.png', dpi=700)

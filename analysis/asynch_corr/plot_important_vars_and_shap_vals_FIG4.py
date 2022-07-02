import numpy as np
import rioxarray as rxr
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib import ticker
import palettable
import os


# PARAMS:
# percentiles (expressed as diff from 0 and 100)
pctile_offset = 1
# colormaps
cmap_asynch_covar = 'inferno'
cmap_other_covar = 'cividis'
covar_cmap_dict = {True: cmap_asynch_covar, False: cmap_other_covar}
cmap_shap = palettable.cmocean.diverging.Curl_20.mpl_colormap
#cmap_shap = palettable.scientific.diverging.Berlin_20.mpl_colormap
# fontsizes
title_fontsize = 26
ylab_fontsize = 22
cbar_ticklab_fontsize = 14
# save fig?
save_plot = True

# set directories
data_dir = '/media/deth/SLAB/seasonality/'
plot_dir = '/home/deth/Desktop/CAL/research/projects/seasonality/results/figs'

# NOTE: somewhat arbitrary, but just showing 3 highest-import variables

# filenames for SHAP maps and original data
vars = {'ppt.asy': 'pr_global_asynch.tif',
        'ppt.sea.nsd': 'CHELSA_bio15_1981-2010_V.2.1_5km_10CELLRAD_NEIGHSD.tif',
        'tmp.min.asy': 'tmmn_global_asynch.tif',
        #'tmp.min.nsd': 'CHELSA_bio6_1981-2010_V.2.1_5km_10CELLRAD_NEIGHSD.tif',
       }

# get country boundaries
countries = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
# reproject to equal-earth
countries = countries.to_crs(8857)


# function for formatting axes
def format_axes(ax, ylabel, bounds, im):
    # add countries
    countries.plot(color='none',
                   edgecolor='black',
                   linewidth=0.5,
                   zorder=1,
                   ax=ax)
    # clear axis text and labels (except setting y-labels on left side)
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_xlabel('')
    ax.set_ylabel(ylabel,
                  fontdict={'fontsize': ylab_fontsize})
    ax.set_title('')
    # set axis bounds
    xlim = bounds[::2]
    ylim = bounds[1::2]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    # NOTE: CLIP OUT HAWAII AND MUCH OF PACIFIC, TO BETTER FOCUS MAP,
    #       AND RE-INCLUDE SOUTHERN TIP OF S. AMERICA
    ax.set_xlim([ax.get_xlim()[0]*0.84, ax.get_xlim()[1]])
    ax.set_ylim([ax.get_ylim()[0]*1.06, ax.get_ylim()[1]])
    # add and format colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=13, length=0) # labelsize, no ticks
    tick_locator = ticker.MaxNLocator(nbins=5)  # 5 ticks
    cbar.locator = tick_locator
    if ylabel != '':
        cbar.formatter.set_powerlimits((0, 0)) # scientific notation on left col
    cbar.update_ticks()

    return cbar


# set up figure
fig = plt.figure(figsize=(15, 3.3*len(vars)))
gs = fig.add_gridspec(len(vars),2)

# set row counter
row = 0

# plot each var's maps
for var, covar_file in vars.items():
    print('\nplotting %s...\n\n' % var)

    # load and reproject Shapley-value map
    ax_shap = fig.add_subplot(gs[row, 1])
    shap_file = 'SHAP_map_' + var + '.tif'
    rast_shap = rxr.open_rasterio(os.path.join(data_dir,
                                          'results',
                                          shap_file), masked=True)[0]
    rast_shap = rast_shap.rio.reproject(countries.crs)
    # get rid of anomalous huge values
    rast_shap.values[rast_shap.values > 1e300] = np.nan
    # mask to match footprint of country boundaries
    #rast_shap = rast_shap.rio.clip(countries.geometry)


    # load and reproject covariate map
    ax_covar = fig.add_subplot(gs[row, 0])
    rast_covar = rxr.open_rasterio(os.path.join(data_dir,
                                          'other/rf_vars',
                                          covar_file), masked=True)
    # NOTE: need to grab the right layer number
    if 'asy' in var:
        rast_covar = rast_covar[2]
    else:
        rast_covar = rast_covar[0]
    rast_covar = rast_covar.rio.reproject(countries.crs)
    # get rid of anomalous huge values
    rast_shap.values[rast_shap.values > 1e300] = np.nan
    # mask to match footprint of country boundaries
    rast_covar = rast_covar.rio.clip(countries.geometry)

    # plot both maps
    try:
        vmin_covar, vmax_covar = np.nanpercentile(rast_covar.values, [pctile_offset,
                                                          100-pctile_offset])
        # set asynchrony viz min value to 0
        # NOTE: I checked this against pctile plots and makes plenty of sense
        if 'asy' in var:
            vmin_covar = 0
        im_covar = rast_covar.plot.imshow(ax=ax_covar,
                                          #cmap = covar_cmap_dict['asy' in var],
                                          cmap=cmap_other_covar,
                                          add_colorbar=False,
                                          zorder=0,
                                          vmin=vmin_covar,
                                          vmax=vmax_covar,
                                         )
        # NOTE: JUST FIXING ALL SHAP MAPS AT THE SAME MIN-MAX VALUES
        #       (DRAWN EMPIRICALLY FROM MAP OF TOP-IMPORTANCE VARIABLE: ppt.asy)
        #       TO GIVE RELATIVE IDEA OF IMPORTANCE OF PREDICTED VALS 
        vmin_shap, vmax_shap = (-0.65, 0.65)
        im_shap = rast_shap.plot.imshow(ax=ax_shap,
                                        cmap=cmap_shap,
                                        add_colorbar=False,
                                        zorder=0,
                                        vmin=vmin_shap,
                                        vmax=vmax_shap,
                                       )
    # NOTE: except AttributeError from broken handling of potential LaTeX text
    except AttributeError:
        pass

    # format both plots' axes
    bounds = rast_shap.rio.bounds()
    cbar_covar = format_axes(ax_covar, var, bounds, im_covar)
    cbar_shap = format_axes(ax_shap, '', bounds, im_shap)
    # redo SHAP map cbar ticks
    cbar_shap.ax.set_yticks([vmin_shap, vmin_shap/2, 0, vmax_shap/2, vmax_shap])

    # add col titles at top of first row
    if row == 0:
        ax_covar.set_title('covariate',
                           fontdict={'fontsize': title_fontsize})
        ax_shap.set_title('SHAP values',
                          fontdict={'fontsize': title_fontsize})

    #increment row counter
    row += 1

# adjust subplots, show, and save
plt.subplots_adjust(left=0.06,
                    right=0.92,
                    bottom=0.02,
                    top=0.94,
                    hspace=0.1,
                    wspace=0.1)
plt.show()
if save_plot:
    fig.savefig(os.path.join(plot_dir, 'fig4_rf_results.png'), dpi=700)

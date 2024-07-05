import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
from rasterstats import zonal_stats
from rasterio.enums import Resampling
import os
import sys
import re


sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf

# set data-directory path
data_dir = phf.EXTERNAL_MASK_DATA_DIR

# general plotting params:
title_fontsize = 12
rowlab_fontsize = 18
axislab_fontsize = 13
ticklab_fontsize = 12
cbarlab_fontsize = 12
cbar_axlab_fontsize = 35
cbar_ticklab_fontsize = 24
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
supp_axlabel_fontdict = {'fontsize': 18}
cbarlabel_fontdict = {'fontsize': 14}
cbarticklabel_fontsize = 12
partlabel_fontsize=40
target_crs = 8857


# load the water mask (just for nicer plotting)
water_mask_unproj = rxr.open_rasterio(os.path.join(data_dir,
                                            'waterMask.tif'), masked=False)[0]
water_mask= water_mask_unproj.rio.write_crs(4326).rio.reproject(target_crs,
                                                  resampling=Resampling.mode,
                                                 )

# load world continents dataset from
# https://github.com/Esri/arcgis-runtime-samples-data/tree/master/shapefiles
# to use for map tabulation by continent
cont = gpd.read_file(os.path.join(phf.BOUNDS_DIR, 'world-continents.shp'))
cont = cont[cont['CONTINENT'] != 'Antarctica']
new_oceania = cont[cont['CONTINENT'].isin(['Australia',
                                           'Oceania'])].dissolve()['geometry']
cont.loc[cont['CONTINENT'] == 'Oceania', 'geometry'] = new_oceania.values[0]
cont = cont[cont['CONTINENT'] != 'Australia']


def map_mask(ax, mask_filename, axlabel, lcMask_mode=None):
    """
    plot a masking map from the LSP-fitting process
    """
    files = [f for f in os.listdir(data_dir) if
             os.path.split(f)[-1] == mask_filename]
    assert len(files) == 1
    mask = rxr.open_rasterio(os.path.join(data_dir, files[0]),
                             masked=True)[0]
    # NOTE: for some reason, the short ts mask and significance mask
    #       saved with one extra row and column and with order of y coordinates
    #       reversed, but a simple reproject_match to the unprojected
    #       water_mask fixes this
    assert mask.shape in [(2700, 6900), (2701, 6901)]
    if mask.shape == (2701, 6901):
        mask = mask.rio.reproject_match(water_mask_unproj)
    if 'lcMask' in mask_filename:
        unique_vals = np.array((0, 1, 2))
    else:
        unique_vals = np.array((0, 1))
    assert np.all(np.unique(mask.values) == unique_vals)
    # get zonal stats
    if 'lcMask' in mask_filename:
        # get zonal stats for both 'default' and 'strict' land cover masking
        # (i.e., mask for all values >0 (default) and all values >1 (strict))
        zonal_means = [zonal_stats(vectors=cont.to_crs(mask.rio.crs),
                                   raster=(mask > val).values,
                                   affine=mask.rio.transform(),
                                   stats=["mean"],
                                   geojson_out=True,
                                  ) for val in range(2)]
    else:
        zonal_means = [zonal_stats(vectors=cont.to_crs(mask.rio.crs),
                                   raster=mask.values,
                                   affine=mask.rio.transform(),
                                   stats=["mean"],
                                   geojson_out=True,
                                  )]
    mask_proj = mask.rio.write_crs(4326).rio.reproject(target_crs,
                                                  resampling=Resampling.mode,
                                                 )
    # get rid of strips of large values from reprojection
    mask_proj = mask_proj.where(mask_proj<=2, np.nan)
    # NOTE: annoying AttributeError is because da.attrs['long_name']
    #       is retained as a tuple of names (rather than being subsetted
    #       by indexing) when I index a single layer out of an
    #       xarray.core.dataarray.DataArray;
    #       for now, a hacky fix is just assigning a string to that attr
    mask_proj.attrs['long_name'] = ''
    if 'lcMask' in mask_filename:
        colors = ['red', 'black', 'white']
        vmax = 2
    else:
        colors = ['red', 'white']
        vmax = 1
    cmap = ListedColormap(colors)
    # mask out oceans, seas, and other large water bodies
    mask_proj = mask_proj.where((water_mask==1).values)
    # mask outside bounds of Equal Area projection
    mask_proj = phf.mask_xarr_to_other_xarr_bbox(mask_proj, mask)
    mask_proj.plot.imshow(ax=ax,
                          cmap=cmap,
                          vmin=0,
                          vmax=vmax,
                          add_colorbar=False,
                          zorder=0,
                         )
    phf.plot_juris_bounds(ax,
                          lev0_color='#060606',
                          lev0_linecolor='gray',
                          lev0_linewidth=0.5,
                          lev0_alpha=0.8,
                          lev0_zorder=2,
                          lev1_color='#060606',
                          lev1_linecolor='gray',
                          lev1_linewidth=0.3,
                          lev1_alpha=0.7,
                          lev1_zorder=1,
                          strip_axes=True,
                          crs=target_crs,
                         )
    ax.set_xlim(mask_proj.rio.bounds()[0::2])
    ax.set_ylim(mask_proj.rio.bounds()[1::2])
    ax.text(-1.4e7, -6.4e6, axlabel, size=7)
    # process the zonal means and return them
    # NOTE: taking 1 minus the calculated mean, to express as percent masked out
    zonal_means = [{j['properties']['CONTINENT']: 1-j['properties']['mean']
                    for j in geojson} for geojson in zonal_means]
    return zonal_means


if __name__ == '__main__':

    # dict to store masking fractions for the summary table
    fracs = {c: [] for c in cont['CONTINENT'].values}
    fracs['mask'] = []

    # dict for renaming masks, for display in fig and in table
    label_dict = {'lcMask': 'land cover',
                  'monthPropsMinMask_NIRv': 'monthly data avail.',
                  'evennessMask_NIRv': 'monthly data evenness',
                  'shortTSMask_NIRv': 'total data avail.',
                  'signifMask_NIRv': 'regression signif.',
                 }

    fig = plt.figure(figsize=(8.7, 6))
    gs = GridSpec(3, 2, figure=fig)
    filenames = ['lcMask.tif',
                 'monthPropsMinMask_NIRv.tif',
                 'evennessMask_NIRv.tif',
                 'shortTSMask_NIRv.tif',
                 'signifMask_NIRv.tif',
                ]

    for ct, fn in enumerate(filenames):
        print(f"\n\nNOW PROCESSING {fn}...\n\n")
        ax = fig.add_subplot(gs[ct//2, ct%2])
        label = label_dict[fn.replace('.tif', '')]
        zonal_means = map_mask(ax, fn, label)
        if 'lcMask' in fn:
            fracs['mask'].append(label+' (default mask)')
            fracs['mask'].append(label+' (strict mask)')
        else:
            fracs['mask'].append(label)
        [fracs[c].append(v) for zm in zonal_means for c, v in zm.items()]

    fig.subplots_adjust(left=0,
                        right=1,
                        bottom=0.01,
                        top=0.99,
                        hspace=0,
                        wspace=0,
                       )
    fig.savefig(os.path.join(phf.FIGS_DIR, 'FIG_SUPP_mask_maps.png'), dpi=500)

    # save the mask-fraction table
    frac_tab = pd.DataFrame.from_dict(fracs)
    frac_tab = frac_tab.iloc[:, [6]+[*range(6)]]
    frac_tab.to_csv(os.path.join(phf.TABS_DIR,
                                 'TAB_SUPP_mask_fracs.csv'),
                    index=False,
                   )


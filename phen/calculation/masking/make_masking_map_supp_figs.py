import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
from rasterstats import zonal_stats
from rasterio.enums import Resampling
import os
import sys
import re


# TODO:
    # - GET RID OF LINE CHANGING SINGLE PIXEL IN SUMMED LC MASK
    # - DELETE '_DUMMY' IN LABEL-PROCESSING CODE



sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf

# set data-directory path
data_dir = phf.EXTERNAL_RF_DATA_DIR

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

# load shapefiles
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                    'NewWorldFile_2020.shp')).to_crs(target_crs)
# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                                  f)).to_crs(target_crs))
subnational = pd.concat(subnational)

# load world continents dataset from
# https://github.com/Esri/arcgis-runtime-samples-data/tree/master/shapefiles
# to use for map tabulation by continent
cont = gpd.read_file(os.path.join(phf.BOUNDS_DIR, 'world-continents.shp'))
cont = cont[cont['CONTINENT'] != 'Antarctica']
new_oceania = cont[cont['CONTINENT'].isin(['Australia',
                                           'Oceania'])].dissolve()['geometry']
cont.loc[cont['CONTINENT'] == 'Oceania', 'geometry'] = new_oceania.values[0]
cont = cont[cont['CONTINENT'] != 'Australia']

def plot_juris_bounds(ax, subnat_zorder=0, nat_zorder=1,
                      polys_color='#060606',
                     ):
    """
    plot national and subnational jurisdictional bounds
    """
    subnational.plot(ax=ax,
                     color='none',
                     edgecolor='gray',
                     linewidth=0.3,
                     alpha=0.7,
                     zorder=subnat_zorder,
                    )
    countries.plot(ax=ax,
                   color='none',
                   edgecolor='gray',
                   linewidth=0.5,
                   alpha=0.8,
                   zorder=nat_zorder,
                  )


def map_mask(ax, mask_filename, axlabel):
    """
    plot a masking map from the LSP-fitting process
    """
    if isinstance(mask_filename, str):
        files = [f for f in os.listdir(phf.EXTERNAL_DATA_DIR) if
                 os.path.split(f)[-1] == mask_filename]
        assert len(files) == 1
        mask = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR, files[0]),
                             masked=True)[0]
        assert np.all(np.unique(mask.values) == np.array((0, 1)))
        # get zonal stats
        zonal_means = [zonal_stats(vectors=cont.to_crs(mask.rio.crs),
                                   raster=mask.values,
                                   affine=mask.rio.transform(),
                                   stats=["mean"],
                                   geojson_out=True,
                                  )]
    elif isinstance(mask_filename, list):
        files = [f for f in os.listdir(phf.EXTERNAL_DATA_DIR) if
                 os.path.split(f)[-1] in mask_filename]
        assert len(files) == 2
        mask = [rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR, f),
                                  masked=True)[0] for f in files]
        zonal_means = [zonal_stats(vectors=cont.to_crs(m.rio.crs),
                                   raster=m.values,
                                   affine=m.rio.transform(),
                                   stats=["mean"],
                                   geojson_out=True,
                                  ) for m in mask]
        # TODO: DELETE ME
        mask[0][10:100, 10:100] = 1
        mask = mask[0] + mask[1]
        assert np.all(np.unique(mask.values) == np.array((0, 1, 2)))
    mask = mask.rio.write_crs(4326).rio.reproject(target_crs,
                                                  resampling=Resampling.mode,
                                                 )
    # get rid of strips of large values from reprojection
    mask = mask.where(mask<=2, np.nan)
    # NOTE: annoying AttributeError is because da.attrs['long_name']
    #       is retained as a tuple of names (rather than being subsetted
    #       by indexing) when I index a single layer out of an
    #       xarray.core.dataarray.DataArray;
    #       for now, a hacky fix is just assigning a string to that attr
    mask.attrs['long_name'] = ''
    if isinstance(mask_filename, str):
        cmap = 'gist_gray'
        vmax = 1
    elif isinstance(mask_filename, list):
        colors = ['black', 'red', 'white']
        cmap = ListedColormap(colors)
        vmax = 2
    mask.plot.imshow(ax=ax,
                     cmap=cmap,
                     vmin=0,
                     vmax=vmax,
                     add_colorbar=False,
                     zorder=0,
                    )
    plot_juris_bounds(ax=ax)
    ax.set_xlim(mask.rio.bounds()[0::2])
    ax.set_ylim(mask.rio.bounds()[1::2])
    # NOTE: chopping off western edge because the equal earth projection
    #       makes NZ appear twice
    ax.set_xlim(0.95 * ax.get_xlim()[0], ax.get_xlim()[1])
    ax.set_title(axlabel, fontdict={'fontsize': 10})
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(())
    ax.set_yticks(())
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
    label_dict = {'lcSummaryMask_NIRv_default': 'land cover',
                  'lcSummaryMask_NIRv_strict': 'land cover',
                  'lcStabilityMask_NIRv': 'land cover stability',
                  'monthPropsMinMask_NIRv': 'monthly data avail.',
                  'evennessMask_NIRv': 'monthly data evenness',
                  'shortTSMask_NIRv': 'total data avail.',
                  'signifMask_NIRv': 'regression signif.',
                 }

    fig = plt.figure(figsize=(8, 11))
    filenames = [['lcSummaryMask_NIRv_default_DUMMY.tif',
                  'lcSummaryMask_NIRv_strict_DUMMY.tif'],
                 'lcStabilityMask_NIRv_DUMMY.tif',
                 'monthPropsMinMask_NIRv_DUMMY.tif',
                 'evennessMask_NIRv_DUMMY.tif',
                 'shortTSMask_NIRv_DUMMY.tif',
                 'signifMask_NIRv_DUMMY.tif',
                ]

    for ct, fn in enumerate(filenames, 1):
        print(f"\n\nNOW PROCESSING {fn}...\n\n")
        ax = fig.add_subplot(3,2,ct)
        if isinstance(fn, list):
            # TODO: DELETE '_DUMMY' BELOW
            label = label_dict[fn[0].replace('_DUMMY.tif', '')]
            zonal_means = map_mask(ax, fn, label)
            for i, f in enumerate(fn):
                if 'strict' in f:
                    sublabel = label + ' (incl. ag.)'
                else:
                    sublabel = label
                fracs['mask'].append(sublabel)
                [fracs[c].append(v) for c, v in zonal_means[i].items()]
        elif isinstance(fn, str):
            # TODO: DELETE '_DUMMY' BELOW
            label = label_dict[fn.replace('_DUMMY.tif', '')]
            zonal_means = map_mask(ax, fn, label)
            fracs['mask'].append(label)
            [fracs[c].append(v) for c, v in zonal_means[0].items()]


    fig.subplots_adjust(left=0.02,
                        right = 0.98,
                        bottom=0.08,
                        top=0.98,
                        hspace=0.8,
                        wspace=0.0,
                       )
    fig.savefig('FIG_SXXX_mask_maps.png', dpi=500)

    # save the mask-fraction table
    frac_tab = pd.DataFrame.from_dict(fracs)
    frac_tab = frac_tab.iloc[:, [6]+[*range(6)]]
    frac_tab.to_csv('TAB_SXXX_mask_fracs.csv', index=False)


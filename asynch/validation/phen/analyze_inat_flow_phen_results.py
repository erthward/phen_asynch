#!/usr/bin/env python
# coding: utf-8


##########################################################################
# TODO:

    # what to do about species with no clear seasonality signal? permut test?

    # what to do about species with huge range across N and S hemispheres, or
    # across numerous continents?

    # figure out way forward for environmental regssions

    # develop better way of determining/categorizing PGF

    # consider using GIFT database to correct for proportion of regional taxa
    # assessed

##########################################################################


import os
import re
import h3
import sys
import rasterstats
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely.geometry import Polygon, Point, MultiPolygon


#------------------------------------------------------------------------------
# params:
#------------------------------------------------------------------------------

# minimum significant R^2ratio
min_r2ratio = 2

# max dip p-value
max_dip_pval = 0.05

# directory where hex GeoJSON is to be saved
hex_data_dir = '/media/deth/SLAB/diss/3-phn/inat/'

# directory with downloaded iNat observation datasets
obs_data_dir = os.path.join(hex_data_dir, 'obs_data')


#------------------------------------------------------------------------------
# load data:
#------------------------------------------------------------------------------

# load iNat phenology results
res_gdf = gpd.read_file('inat_flower_phen_results.json')
len_b4 = len(res_gdf)
res_gdf = res_gdf[(pd.notnull(res_gdf['r2_ratio']) &
                   pd.notnull(res_gdf['dip_pval']) &
                   pd.notnull(res_gdf['geometry']))]
len_af = len(res_gdf)
print((f"{len_b4-len_af} rows dropped because of missing geometries or "
        "statistical results."))

# merge on TRY plant growth forms
all_tax_file = './all_inat_plant_phen_taxa_w_TRY_pgf.csv'
cols = ['tid', 'pgf']
res_gdf = pd.merge(res_gdf,
                   pd.read_csv(all_tax_file).loc[:, cols],
                   on='tid',
                   how='left',
                  )
res_gdf['tree'] = [(re.search('w', s) is not None and
                       re.search('f', s) is None) if pd.notnull(s) else False
                                                for s in res_gdf['pgf'].values]
res_gdf['nontree'] = [re.search('w', s) is None if pd.notnull(s) else False
                                                for s in res_gdf['pgf'].values]

# local modules
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf

# set common equal-area projection to use
crs = 8857

# load country boundaries
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                       'NewWorldFile_2020.shp')).to_crs(crs)

# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR,f)).to_crs(crs))
subnational = pd.concat(subnational)

# load asynch data (band 0 in the asynch file
neigh_rad = 100
asynch = rxr.open_rasterio(phf.ASYNCH_FILES[neigh_rad])[0]

# load climate data
bioclim = [rxr.open_rasterio(f, masked=True) for f in phf.BIOCLIM_INFILEPATHS]


#------------------------------------------------------------------------------
# create hexes:
#------------------------------------------------------------------------------
hex_filename = 'inat_hex.json'
if not os.path.isfile(os.path.join(hex_data_dir, hex_filename)):
    # NOTE: based on code taken from:
        # https://towardsdatascience.com/uber-h3-for-data-analysis-with-python-1e54acdcc908
    # make dataframe to hold h3-converted data
    h3_df = pd.DataFrame([], columns=['row_id', 'h3_id',
                                      'h3_geo_boundary', 'h3_centroid'])
    # loop over results rows and convert to H3 hexes
    for i, row in res_gdf.iterrows():
        p = row['geometry']
        if isinstance(p, MultiPolygon):
            ps = list(p.geoms)
        else:
            ps = [p]
        for poly in ps:
            poly_json = gpd.GeoSeries([poly]).__geo_interface__
            poly_json = poly_json['features'][0]['geometry']
            h3_hexes = h3.polyfill_geojson(poly_json, 3)
            for h3_hex in h3_hexes:
                h3_geo_boundary = Polygon(
                    h3.h3_to_geo_boundary(h3_hex,geo_json=True))
                h3_centroid = h3.h3_to_geo(h3_hex)
                h3_df.loc[len(h3_df)]= [i, h3_hex, h3_geo_boundary, h3_centroid]
    # coerce to GeoDataFrame
    geoms = [Polygon(row['h3_geo_boundary']) for i, row in h3_df.iterrows()]
    h3_df['geometry'] = geoms
    h3_gdf = gpd.GeoDataFrame(h3_df, geometry='geometry', crs=4326)
    # deduplicate hexes
    h3_gdf = h3_gdf.drop_duplicates(subset='geometry')
    # rework columns so that we can save to GeoJSON
    h3_gdf['x'] = [c[1] for c in h3_gdf['h3_centroid'].values]
    h3_gdf['y'] = [c[0] for c in h3_gdf['h3_centroid'].values]
    h3_gdf = h3_gdf.loc[:, ['row_id', 'h3_id', 'x', 'y', 'geometry']]
    h3_gdf.to_file(os.path.join(hex_data_dir, hex_filename))
else:
    h3_gdf = gpd.read_file(os.path.join(hex_data_dir, hex_filename))


#------------------------------------------------------------------------------
# summarize iNat phen results within hexes:
#------------------------------------------------------------------------------
res_hex_filename = 'inat_hex_results.json'
if not os.path.isfile(os.path.join(hex_data_dir, res_hex_filename)):
    # summarize results within hexes
    new_cols = {'n_taxa': [],
                'mean_r2ratio': [],
                'mean_dip_pval': [],
                'prop_r2ratio_signif': [],
                'prop_dip_pval_signif': [],
                'prop_both_signif': [],
                'tree_n_taxa': [],
                'tree_prop_r2ratio_signif': [],
                'tree_prop_dip_pval_signif': [],
                'tree_prop_both_signif': [],
                'nontree_n_taxa': [],
                'nontree_prop_r2ratio_signif': [],
                'nontree_prop_dip_pval_signif': [],
                'nontree_prop_both_signif': [],
               }
    tree_gdf = res_gdf[res_gdf['tree']]
    nontree_gdf = res_gdf[res_gdf['nontree']]
    print(f"\n\n{len(tree_gdf)} tree-annotated taxa\n\n")
    print(f"\n\n{len(nontree_gdf)} nontree-annotated taxa\n\n")
    pgf_gdfs = {'tree': tree_gdf, 'nontree': nontree_gdf}
    for i, row in h3_gdf.iterrows():
        # get all taxa that intersect with this hex
        gdf_intsxn = res_gdf[res_gdf.intersects(row['geometry'])]
        new_cols['n_taxa'].append(len(gdf_intsxn))
        new_cols['mean_r2ratio'].append(np.mean(gdf_intsxn['r2_ratio']))
        new_cols['mean_dip_pval'].append(np.mean(gdf_intsxn['dip_pval']))
        r2ratio_signif = gdf_intsxn['r2_ratio'] >= min_r2ratio
        dip_pval_signif = gdf_intsxn['dip_pval'] <= max_dip_pval
        both_signif = (r2ratio_signif) & (dip_pval_signif)
        new_cols['prop_r2ratio_signif'].append(np.mean(r2ratio_signif))
        new_cols['prop_dip_pval_signif'].append(np.mean(dip_pval_signif))
        new_cols['prop_both_signif'].append(np.mean(both_signif))
        # do the same for only trees and only non-trees
        for pgf, pgf_gdf in pgf_gdfs.items():
            pgf_gdf_intsxn = pgf_gdf[pgf_gdf.intersects(row['geometry'])]
            new_cols[f'{pgf}_n_taxa'].append(len(pgf_gdf_intsxn))
            r2ratio_signif = pgf_gdf_intsxn['r2_ratio'] >= min_r2ratio
            dip_pval_signif = pgf_gdf_intsxn['dip_pval'] <= max_dip_pval
            both_signif = (r2ratio_signif) & (dip_pval_signif)
            new_cols[f'{pgf}_prop_r2ratio_signif'].append(np.mean(r2ratio_signif))
            new_cols[f'{pgf}_prop_dip_pval_signif'].append(np.mean(dip_pval_signif))
            new_cols[f'{pgf}_prop_both_signif'].append(np.mean(both_signif))
    for col, vals in new_cols.items():
        h3_gdf[col] = vals
    for c in new_cols.keys():
        if c.startswith('prop'):
            assert np.nanmin(h3_gdf[c]) >= 0 and np.nanmax(h3_gdf[c]) <= 1
    # merge asynch values onto the hex gdf
    assert asynch.rio.crs.wkt == h3_gdf.crs
    zon_stats = rasterstats.zonal_stats(vectors=h3_gdf['geometry'],
                                        raster=asynch.values,
                                        affine=asynch.rio.transform(),
                                        nodata=-9999,
                                        stats=['mean', 'median'],
                                       )
    h3_gdf['mean_asynch'] = [zs['mean'] for zs in zon_stats]
    h3_gdf['medn_asynch'] = [zs['median'] for zs in zon_stats]
    h3_gdf.to_file(os.path.join(hex_data_dir, res_hex_filename))
else:
    h3_gdf = gpd.read_file(os.path.join(hex_data_dir, res_hex_filename))


#------------------------------------------------------------------------------
# plot results
#------------------------------------------------------------------------------
fig = plt.figure(figsize=(9,15))
cmap='viridis'
# plot results for trees, non-trees, and all taxa combined
res_cols = ['tree_prop_both_signif',
            'nontree_prop_both_signif',
            'prop_both_signif',
           ]
label_dict = {0: 'trees',
              1: 'non-trees',
              2: 'all taxa',
             }
for i, res_col in enumerate(res_cols):
    ax = fig.add_subplot(3, 1, i+1)
    # add bottom axes for a colorbar
    if i == 2:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('bottom', size='7%', pad=0.2)
    # transform to equal-area projection and plot
    subnational.plot(color='none',
                     edgecolor='black',
                     zorder=0,
                     ax=ax,
                     alpha=0.6,
                     linewidth=0.05,
                    )
    countries.plot(color='none',
                                edgecolor='black',
                                linewidth=0.1,
                                zorder=1,
                                ax=ax,
                                )
    h3_gdf.to_crs(crs).plot(res_col,
                            cmap=cmap,
                            alpha=0.75,
                            zorder=2,
                            ax=ax,
                            edgecolor='white',
                            linewidth=0.2,
                            legend=False,
                            vmin=0,
                            #vmax=1,
                            vmax=np.percentile(h3_gdf['prop_both_signif'], 95),
                           )
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_title('')
    # clip off longitudinally and latitudinally
    #ax.set_xlim(0.85 * ax.get_xlim()[0], 0.87*ax.get_xlim()[1])
    #ax.set_ylim(1.11*np.min(h3_gdf.to_crs(crs).geometry.centroid.y),
    #            1.175*np.max(h3_gdf.to_crs(crs).geometry.centroid.y))
    if i  == 2:
        scalcmap = plt.cm.ScalarMappable(cmap=cmap,
                                         norm=plt.Normalize(vmin=0, vmax=1),
                                        )
        plt.colorbar(scalcmap, cax=cax, orientation='horizontal')
        xticks = np.linspace(0, 1, 5)
        cax.set_xlabel('proportion significantly non-unimodal taxa',
                       fontdict={'fontsize': 11},
                      )
        cax.set_xticks(xticks, ['%0.2f' % t for t in xticks], size=9)
        cax.set_ylabel('')
        cax.set_yticks(())
# adjust subplots and write to disk
fig.subplots_adjust(top=0.96,
                    bottom=0.125,
                    right=0.96,
                    left=0.04,
                    wspace=0.5,
                    hspace=0.3,
                   )
#fig.savefig(...

################################
################################
################################
################################
################################
################################

assert False


fig, ax = plt.subplots(1,1)
#subnational.plot(color='none',
#                     edgecolor='black',
#                     zorder=0,
#                     ax=ax,
#                     alpha=0.6,
#                     linewidth=0.05,
#                    )
countries.plot(color='none',
                                edgecolor='black',
                                linewidth=0.1,
                                zorder=1,
                                ax=ax,
                                )
h3_gdf.to_crs(crs).plot('prop_dip_pval_signif',
                            cmap=cmap,
                            alpha=0.75,
                            zorder=2,
                            ax=ax,
                            edgecolor='white',
                            linewidth=0.2,
                            legend=True,
                            vmin=0,
                            vmax=1,
                            #vmax=np.percentile(h3_gdf['prop_both_signif'], 95),
                           )







# focal regressions for certain species
#tid = 58040 # Stephanomeria pauciflora
tid = 327909 # Turnera subulata
files = [f for f in os.listdir(obs_data_dir) if re.search(f"TID_{tid}", f)]
assert len(files) == 1
file = files[0]
obs = gpd.read_file(os.path.join(obs_data_dir, file))
obs['x'] = obs['geometry'].x
obs['y'] = obs['geometry'].y
for bc_i, bc in enumerate(bioclim):
    obs[f'bio_{bc_i}'] = [bc.sel(x=row['x'], y=row['y'],
                    method='nearest').values[0] for i, row in obs.iterrows()]
bc_data = obs.loc[:, [c for c in obs.columns if c.startswith('bio_')]]
complete_rows = np.all(pd.notnull(bc_data), axis=1)
bc_data = bc_data.loc[complete_rows, :]
pcs = PCA(n_components=len(bioclim)).fit(bc_data)
plt.plot(pcs.explained_variance_ratio_)
top_n = 3
pc_vals = pcs.fit_transform(bc_data)[:, :top_n]
reg = LinearRegression().fit(pc_vals, np.sin(obs.loc[complete_rows, :]['doy_circ']))
reg.score(pc_vals, np.sin(obs.loc[complete_rows, :]['doy_circ']))
reg.coef_
plt.scatter(pc_vals[:, 0], np.sin(obs.loc[complete_rows, :]['doy_circ']))



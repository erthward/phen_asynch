#!/usr/bin/env python
# coding: utf-8


##########################################################################
# TODO:

    # biggest reason MMRRs are taking forever, I believe, is the raster
    # interpolation. need to decide:
        # A. do we do it? if so how much?
        # B. if we do, probably much faster to do once, save the interpolated
        #    raster, then read from that, right?

    # what to do about species with huge range across N and S hemispheres, or
    # across numerous continents? drop completely? or not a problem if distance
    # partialled out? or try to cluster regions and analyze separately?

    # develop better way of determining/categorizing PGF

    # consider using GIFT database to correct for proportion of regional taxa
    # assessed

##########################################################################


import os
import re
import h3
import sys
import pyproj
import rasterstats
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import OrderedDict
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely.geometry import Polygon, Point, MultiPolygon

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf
from MMRR import MMRR



#------------------------------------------------------------------------------
# params:
#------------------------------------------------------------------------------

# whether to make maps and whether to run MMRR models
make_maps = False
run_MMRRs = True

# minimum significant R^2ratio
min_r2ratio = 2

# max dip p-value
max_dip_pval = 0.05

# whether to interpolate missing LSP coefficients, 
# and if so, how far away can rasterio look for neighbor to interpolate from?
# coefficient raster?
interp_sea_data = False
neigh_dist_sea_fill_tol = 2

# how many MMRR permutations to use
MMRR_nperm = 999

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
print((f"\n\n{len_b4-len_af} rows dropped because of missing geometries or "
        "missing statistical results."))

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

# calculate a column of our best guess at the number of histogram peaks,
# based on the KDE npeaks and its P-value
res_gdf.loc[:, 'signif_npeaks'] = [row['npeaks'] if row['npeaks_pval']<=0.01
                                   else 0 for i, row in res_gdf.iterrows()]

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
    print("\n\tproducing hex file...\n")
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
    print("\n\treading hexes from file...\n")
    h3_gdf = gpd.read_file(os.path.join(hex_data_dir, hex_filename))


#------------------------------------------------------------------------------
# summarize iNat phen results within hexes:
#------------------------------------------------------------------------------
res_hex_filename = 'inat_hex_results.json'
if not os.path.isfile(os.path.join(hex_data_dir, res_hex_filename)):
    print("\n\tsummarizing results by hex...\n")
    # summarize results within hexes
    new_cols = {'n_taxa': [],
                'mean_r2ratio': [],
                'mean_dip_pval': [],
                'prop_r2ratio_signif': [],
                'prop_dip_pval_signif': [],
                'prop_both_signif': [],
                'prop_1peak_signif': [],
                'prop_0peak_signif': [],
                'prop_2peak_signif': [],
                'tree_n_taxa': [],
                'tree_prop_r2ratio_signif': [],
                'tree_prop_dip_pval_signif': [],
                'tree_prop_both_signif': [],
                'tree_prop_1peak_signif': [],
                'tree_prop_0peak_signif': [],
                'tree_prop_2peak_signif': [],
                'nontree_n_taxa': [],
                'nontree_prop_r2ratio_signif': [],
                'nontree_prop_dip_pval_signif': [],
                'nontree_prop_both_signif': [],
                'nontree_prop_1peak_signif': [],
                'nontree_prop_0peak_signif': [],
                'nontree_prop_2peak_signif': [],
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
        props = []
        for i in range(3):
            prop = np.mean(gdf_intsxn['signif_npeaks'] == i)
            new_cols[f'prop_{i}peak_signif'].append(prop)
            props.append(prop)
        assert np.allclose(np.sum(props), 1)
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
            props = []
            for i in range(3):
                prop = np.mean(pgf_gdf_intsxn['signif_npeaks'] == i)
                new_cols[f'{pgf}_prop_{i}peak_signif'].append(prop)
                props.append(prop)
            assert np.all(pd.isnull(props)) or np.allclose(np.sum(props), 1)
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
    print("\n\treading hex-summarized results from file...\n")
    h3_gdf = gpd.read_file(os.path.join(hex_data_dir, res_hex_filename))


#------------------------------------------------------------------------------
# plot results
#------------------------------------------------------------------------------
if make_maps:
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
                                alpha=1,
                                zorder=2,
                                ax=ax,
                                edgecolor='white',
                                linewidth=0.2,
                                legend=False,
                                vmin=0,
                                #vmax=1,
                                vmax=np.nanpercentile(h3_gdf[res_col], 95),
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
    fig.show()
    #fig.savefig(...



    fig = plt.figure(figsize=(9,15))
    h3_gdf['prop_non1peak_signif'] = (h3_gdf['prop_0peak_signif'] +
                                      h3_gdf['prop_2peak_signif'])
    cmap='viridis'
    # plot results for trees, non-trees, and all taxa combined
    res_cols = ['prop_non1peak_signif',
                'prop_0peak_signif',
                'prop_2peak_signif',
               ]
    label_dict = {0: 'proportion non-unimodal',
                  1: 'proportion with no clear peaks',
                  2: 'proportion with 2 fitted peaks',
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
                                alpha=1,
                                zorder=2,
                                ax=ax,
                                edgecolor='white',
                                linewidth=0.2,
                                legend=False,
                                vmin=0,
                                #vmax=1,
                                vmax=np.nanpercentile(h3_gdf[res_col], 95),
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
    fig.show()


#------------------------------------------------------------------------------
# regress flowering observation doy distance on NIRv LSP data
#------------------------------------------------------------------------------
if run_MMRRs:
    # analyze flow_dist ~ seas_dist for all taxa with non-unimodal flowering hists
    # NOTE: flowering date observations are emphatically NOT flowering period
    #       observations and are inherently noisy, especially for
    #       opportunistically/perennially flowering taxa; if anything, this should
    #       lead us to overlook taxa for which IBT could actually operating on the
    #       ground, and thus should make our overall analysis quite conservative,
    #       such that taxa with significant relationships would be especially good
    #       candidates for potential IBT

    def calc_doy_diff(doy1, doy2):
        '''
        calculate the distance, in number of days, between 2 numericaly days of year
        '''
        d = sorted([doy1, doy2])
        dist = np.min((d[1]-d[0], d[0] + 365 - d[1]))
        return dist


    # set the numpy.random seed
    np.random.seed(1)

    # run analysis for each species with a non-unimodal hist
    # TODO: implement min number of samples here? (because actual n obs usually
    #       < expected n obs because of low positional_accuracy)
    mmrr_dict = {'tid': [],
                 'name': [],
                 'n': [],
                 'r2': [],
                 'int': [],
                 'int_t': [],
                 'int_p': [],
                 'geo': [],
                 'geo_t': [],
                 'geo_p': [],
                 'lsp': [],
                 'lsp_t': [],
                 'lsp_p': [],
                 'f': [],
                 'f_p': [],
                }
    label_dict = {'R^2': 'r2',
                  'Intercept': 'int',
                  'Intercept(t)': 'int_t',
                  'Intercept(p)': 'int_p',
                  'lsp_dist': 'lsp',
                  'lsp_dist(t)': 'lsp_t',
                  'lsp_dist(p)': 'lsp_p',
                  'geo_dist': 'geo',
                  'geo_dist(t)': 'geo_t',
                  'geo_dist(p)': 'geo_p',
                  'F-statistic': 'f',
                  'F p-value': 'f_p',
                 }
    for i, row in res_gdf[res_gdf['signif_npeaks'] != 1].iterrows():
        try:
            tid = row['tid']
            name = row['name']
            print(f'\n\n{"-"*80}\n')
            print(f"running LSP MMRR on flowering observations for {name} (TID: {tid})...")
            # load observation points
            obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
            obs = gpd.read_file(os.path.join(obs_data_dir, obs_fn))

            # get their LSP time series distances
            print('\n\tcalculating LSP distances...')
            pts = np.vstack([obs.geometry.x, obs.geometry.y]).T
            assert pts.shape == (len(obs), 2)
            # TODO: DECIDE WHETHER TO USE DEFAULT OR STRICT COEFFS FILE
            lsp_dist = phf.get_raster_info_points(phf.COEFFS_FILE,
                                                  pts,
                                                  'ts_pdm',
                                                  standardize=True,
                                                  fill_nans=interp_sea_data,
                                                  fill_tol=neigh_dist_sea_fill_tol,
                                                 )
            print((f"\n\n\t{np.round(100*np.mean(pd.isnull(lsp_dist)), 2)} % "
                    "of lsp_dist consists of missing values"))
            keep_idxs = np.where(np.sum(pd.notnull(lsp_dist), axis=1) > 1)[0]
            assert np.all(keep_idxs ==
                          np.where(np.sum(pd.notnull(lsp_dist), axis=0) > 1)[0])
            lsp_dist = lsp_dist[keep_idxs, :][:, keep_idxs]
            assert np.all(lsp_dist.shape == (np.ones(2)*keep_idxs.size))
            assert np.all(lsp_dist == lsp_dist.T)

            # calculate flowering-date distances
            # (i.e., numbers of days between observed flowering dates)
            print('\n\tcalculating flowering-date distances...')
            flw_dist = np.ones(lsp_dist.shape)*np.nan
            keep_obs = obs.iloc[keep_idxs, :]
            for i in range(len(keep_obs)):
                doy_i = keep_obs.iloc[i, :]['doy']
                for j in range(len(keep_obs)):
                    if i==j:
                        flw_dist[i, j] = 0
                    else:
                        doy_j = keep_obs.iloc[j, :]['doy']
                        dist = calc_doy_diff(doy_i, doy_j)
                        flw_dist[i, j] = flw_dist[j, i] = dist
            assert np.sum(pd.isnull(flw_dist)) == 0
            assert np.all(flw_dist == flw_dist.T)

            # calculate geographic distances
            geo_dist = np.ones(flw_dist.shape) * np.nan
            print('\n\tcalculating geographic distances...')
            g = pyproj.Geod(ellps='WGS84')
            keep_pts = pts[keep_idxs, :]
            for i in range(geo_dist.shape[0]):
                lon_i, lat_i = keep_pts[i, :]
                for j in range(geo_dist.shape[1]):
                    # calculate geographic distance
                    if i == j:
                        dist_ij = 0
                    else:
                        lon_j, lat_j = keep_pts[j, :]
                        (az_ij, az_ji, dist_ij) = g.inv(lon_i, lat_i, lon_j, lat_j)
                    geo_dist[i, j] = geo_dist[j, i] = dist_ij
            assert np.sum(pd.isnull(geo_dist)) == 0
            assert np.all(geo_dist == geo_dist.T)

            # run MMRR and store results
            # run the MMRR model and print results
            print('\n\trunning MMRR model...')
            res = MMRR(Y=lsp_dist,
                       X=[geo_dist, flw_dist],
                       Xnames=['lsp_dist', 'geo_dist', 'flw_dist'],
                       # NOTE: MMRR will standardize lower-triangular distance values, and thus
                       #       returns coefficient values as beta-coefficients
                       standardize=True,
                       intercept=True,
                       nperm=MMRR_nperm,
                      )
        except Exception as e:
            print(f"\n\tERROR THROWN: {e}\n\n\tmoving on...\n")
            res = {k: np.nan for k in label_dict.keys()}

        mmrr_dict['tid'].append(tid)
        mmrr_dict['name'].append(name)
        try:
            mmrr_dict['n'].append(lsp_dist.shape[0])
        except Exception:
            mmrr_dict['n'].append(np.nan)
        for k in res.keys():
            mmrr_dict[label_dict[k]].append(res[k])

    # coerce results to pd.DataFrame
    mmrr_df = pd.DataFrame.from_dict(mmrr_dict)
    # sort by increasing P-values of interest and decreasing R^2 within that
    mmrr_df = mmrr_df.sort_values(by=['lsp_p', 'f_p', 'r2'],
                                  ascending=[True, True, False])
    mmrr_df.to_csv('inat_flow_obs_LSP_MMRR_res.csv', index=False)




# load the EOFs, for the plotting function below
eofs = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
            '%s_global_4_EOFs_sqrt_coswts%s%s.tif') %(dataset,
                                                      normts_file_substr,
                                                      mask_filename_ext))[:3]
for i in range(eofs.shape[0]):
    eofs[i] = (eofs[i]-eofs[i].min())/(eofs[i].max()-eofs[i].min())
eofs.rio.set_crs(4326)
eofs = eofs.rio.reproject(crs)
eofs = eofs.where(eofs<3e38, np.nan)

def plot_ts_vs_flow_dates(tid, name,
                          fill_nans=False,
                          fill_tol=None,
                          mmrr_res=None,
                          map=False,
                         ):
    if fill_nans:
        assert fill_tol is not None and isinstance(fill_tol, int)
    if mmrr_res is not None:
        assert ((isinstance(mmrr_res, dict) or
                 isinstance(mmrr_res, OrderedDict)) and
                'lsp_dist' in mmrr_res.keys() and
                'lsp_dist(p)' in mmrr_res.keys())
    obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
    obs = gpd.read_file(os.path.join(obs_data_dir, obs_fn))
    pts = np.vstack([obs.geometry.x, obs.geometry.y]).T
    lsp_ts = phf.get_raster_info_points(phf.COEFFS_FILE,
                                        pts,
                                        'ts',
                                        standardize=True,
                                        fill_nans=fill_nans,
                                        fill_tol=fill_tol,
                                       )
    missing = np.any(pd.isnull(lsp_ts), axis=1)
    print((f"\n{np.round(100*(np.sum(missing)/lsp_ts.shape[0]))}% of sites "
            "are missing LSP data\n"))
    fig = plt.figure()
    ax = fig.add_subplot(1+(int(map)),1,1)
    doy_scale = (obs['doy'].values - np.min(obs['doy']))/(np.max(obs['doy']) -
                                                          np.min(obs['doy']))
    doy_colors = mpl.colormaps.get_cmap('twilight')(doy_scale)
    for i in range(len(obs)):
        if not missing[i]:
            color = doy_colors[i]
            ax.plot(lsp_ts[i, :], color=color)
    title = f"{name} (TID: {tid})"
    if mmrr_res is not None:
        title = (f"{title}\nMMRR LSP coeff: {np.round(mmrr_res['lsp_dist'], 2)}"
                 f" (P={np.round(mmrr_res['lsp_dist(p)'], 2)})")
    ax.set_title(title, fontdict={'fontsize': 12})
    ax.set_xlabel('fitted LSP day of year', fontdict={'fontsize': 12})
    ax.set_ylabel('standardized LSP value', fontdict={'fontsize': 12})
    ax.tick_params(labelsize=10)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    norm = Normalize(vmin=0, vmax=365)
    cbar = ColorbarBase(cax, cmap='twilight', norm=norm, orientation='vertical')
    cbar.set_label('flowering observation day of year', fontsize=12)
    if map:
        ax = fig.add_subplot(1+(int(map)),1,2)
        subnational.plot(color='none',
                         edgecolor='black',
                         zorder=1,
                         ax=ax,
                         alpha=0.6,
                         linewidth=0.05,
                        )
        countries.plot(color='none',
                                    edgecolor='black',
                                    linewidth=0.1,
                                    zorder=2,
                                    ax=ax,
                                    )
        obs.to_crs(crs).plot('doy',
                             cmap='twilight',
                             alpha=1,
                             zorder=3,
                             ax=ax,
                             vmin=0,
                             vmax=365,
                             )
        eofs.plot.imshow(ax=ax,
                         zorder=0,
                        )
        ax.set_xlim(np.min(obs.to_crs(crs).geometry.x), np.max(obs.to_crs(crs).geometry.x))
        ax.set_ylim(np.min(obs.to_crs(crs).geometry.y), np.max(obs.to_crs(crs).geometry.y))
        ax.set_aspect('equal')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_xticks(())
        ax.set_yticks(())
        ax.set_title('')
    return fig

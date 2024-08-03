#!/usr/bin/env python
# coding: utf-8

import os
import re
import h3
import sys
import time
import pprint
import pyproj
import rasterstats
import numpy as np
import pandas as pd
import seaborn as sns
import rasterio as rio
import geopandas as gpd
import rioxarray as rxr
import cmcrameri.cm as cmc
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from collections import OrderedDict
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from shapely.geometry import Polygon, MultiPolygon
from mpl_toolkits.axes_grid1 import make_axes_locatable

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf
from MMRR import MMRR

start_time = time.time()

#------------------------------------------------------------------------------
# params and supporting data:
#------------------------------------------------------------------------------

# whether to interpolate missing LSP coefficients, 
# and if so, how far away can rasterio look for neighbor to interpolate from?
# coefficient raster?
interp_lsp_data = False
neigh_dist_lsp_fill_tol = 2

# how many MMRR permutations to use
# NOTE: using at least enough permutations to make the smallest possible
#       P-value smaller than 0.05 divided by the number of taxa we test
#       after filtering out huge-range taxa,
#       and thus allows for Bonferroni correction
MMRR_nperm = 14999

# common equal-area projection to use
crs = 8857

# whether to use the strict-masked (i.e., ag-removed) coefficients file for
# extracting fitted LSP patterns at flower observation localities
strict_coeffs = True
if strict_coeffs:
    coeffs_file = phf.COEFFS_STRICT_FILE
else:
    coeffs_file = phf.COEFFS_FILE

# directory where hex GeoJSON is to be saved
hex_data_dir = phf.EXTERNAL_INAT_DATA_DIR

# directory with downloaded iNat observation datasets
obs_data_dir = os.path.join(hex_data_dir, 'obs_data')

# load bioclim nodata value
bioclim_nodata_val = rio.open(phf.BIOCLIM_INFILEPATHS[0]).nodata


#------------------------------------------------------------------------------
# load iNaturalist phenology peak-analysis results:
#------------------------------------------------------------------------------
# load iNat phenology results
res_gdf = gpd.read_file('inat_flower_phen_results.json')
len_b4 = len(res_gdf)
res_gdf = res_gdf[pd.notnull(res_gdf['geometry'])]
len_af = len(res_gdf)
fail_fit_hull_msg = (f"\n\n{len_b4-len_af} taxa dropped because of "
                      "failure to fit observation-range alpha hulls.")

# calculate a column of our best guess at the number of histogram peaks,
# based on the KDE npeaks and its P-value
sig_npeaks = [row['npeaks'] if row['npeaks_pval']<=(0.05/len(res_gdf)) else
              0 for i, row in res_gdf.iterrows()]
res_gdf.loc[:, 'signif_npeaks'] = sig_npeaks


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
            poly_json = gpd.GeoSeries([poly]).__geo_interface__['features'][0]['geometry']
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
    print("\n\n\treading hexes from file...\n")
    h3_gdf = gpd.read_file(os.path.join(hex_data_dir, hex_filename))


#------------------------------------------------------------------------------
# summarize iNat phen peak-analysis results within hexes:
#------------------------------------------------------------------------------
peaks_hex_filename = 'inat_hex_results.json'
if not os.path.isfile(os.path.join(hex_data_dir, peaks_hex_filename)):
    print("\n\tsummarizing results by hex...\n")
    # summarize results within hexes
    new_cols = {'n_taxa': [],
                'prop_1peak_signif': [],
                'prop_0peak_signif': [],
                'prop_2pluspeak_signif': [],
               }
    for i, row in h3_gdf.iterrows():
        # get all taxa that intersect with this hex
        gdf_intsxn = res_gdf[res_gdf.intersects(row['geometry'])]
        new_cols['n_taxa'].append(len(gdf_intsxn))
        props = []
        for i in range(3):
            if i < 2:
                prop = np.mean(gdf_intsxn['signif_npeaks'] == i)
                new_cols[f'prop_{i}peak_signif'].append(prop)
            else:
                prop = np.mean(gdf_intsxn['signif_npeaks'] >= i)
                new_cols[f'prop_{i}pluspeak_signif'].append(prop)
            props.append(prop)
        assert np.allclose(np.sum(props), 1)
    for col, vals in new_cols.items():
        h3_gdf[col] = vals
    for c in new_cols.keys():
        if c.startswith('prop'):
            assert np.nanmin(h3_gdf[c]) >= 0 and np.nanmax(h3_gdf[c]) <= 1
    # write it out
    h3_gdf.to_file(os.path.join(hex_data_dir, peaks_hex_filename))
else:
    print("\n\treading hex-summarized phenology peak results from file...\n")
    h3_gdf = gpd.read_file(os.path.join(hex_data_dir, peaks_hex_filename))


#------------------------------------------------------------------------------
# regress flowering observation doy distance on NIRv LSP data
#------------------------------------------------------------------------------
mmrr_res_filename = 'inat_flower_doy_LSP_MMRR_res.csv'
if not os.path.isfile(mmrr_res_filename):

    def fit_LSP_MMRR(obs,
                     strict_coeffs=strict_coeffs,
                     bioclim_nodata_val=bioclim_nodata_val,
                     standardize_lsp_ts=True,
                     interp_lsp_data=interp_lsp_data,
                     neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                     standardize_MMRR_data=True,
                     fit_MMRR_intercept=True,
                     MMRR_nperm=MMRR_nperm,
                    ):
        '''
        fit an MMRR model of the form:
            doy_flower ~ geo_dist + env_dist + LSP_dist
        for the given observation DataFrame
        '''
        # get their LSP time series distances
        print('\n\tcalculating LSP distances...')
        pts = np.vstack([obs.geometry.x, obs.geometry.y]).T
        assert pts.shape == (len(obs), 2)
        lsp_dist = phf.get_raster_info_points(coeffs_file,
                                              pts,
                                              'ts_pdm',
                                              standardize=standardize_lsp_ts,
                                              fill_nans=interp_lsp_data,
                                              fill_tol=neigh_dist_lsp_fill_tol,
                                             )
        pct_lsp_miss = np.mean(pd.isnull(lsp_dist).sum(
                                axis=1) == (lsp_dist.shape[0]-1))
        pct_lsp_miss_msg = (f"{np.round(100*pct_lsp_miss, 2)}% "
                            "of sites are missing LSP coefficients")
        print(pct_lsp_miss_msg)
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
                    dist = phf.calc_doy_diff(doy_i, doy_j)
                    flw_dist[i, j] = flw_dist[j, i] = dist
        assert np.sum(pd.isnull(flw_dist)) == 0
        assert np.all(flw_dist == flw_dist.T)
        # calculate environmental distances
        print('\n\tcalculating environmental distances...')
        keep_pts = pts[keep_idxs, :]
        env_dist = phf.calc_pw_clim_dist_mat(keep_pts,
                                             nodata_val=bioclim_nodata_val,
                                            )
        assert np.sum(pd.isnull(env_dist)) == 0
        assert np.all(env_dist == env_dist.T)
        assert np.all(env_dist.shape == lsp_dist.shape), (('Environmental'
                'distance matrix shape does not match LSP distance matrix'
                                                           'shape!'))
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
        res = MMRR(Y=flw_dist,
                   X=[geo_dist, env_dist, lsp_dist],
                   Xnames=['geo_dist', 'env_dist', 'lsp_dist'],
                   # NOTE: MMRR will standardize lower-triangular distance values, and thus
                   #       returns coefficient values as beta-coefficients
                   standardize=standardize_MMRR_data,
                   intercept=fit_MMRR_intercept,
                   nperm=MMRR_nperm,
                  )
        # get sample size, to return
        n = lsp_dist.shape[0]
        return res, n, pct_lsp_miss


    # analyze flow_dist ~ seas_dist for all taxa with non-unimodal flowering hists
    # NOTE: flowering date observations are emphatically NOT flowering period
    #       observations and thus are inherently noisy, especially for
    #       opportunistically/perennially flowering taxa; if anything, this should
    #       lead us to overlook taxa for which IBT could actually operating on the
    #       ground, and thus should make our overall analysis quite conservative,
    #       such that taxa with significant relationships would be especially good
    #       candidates for potential IBT

    # set the numpy.random seed
    np.random.seed(1)

    # create logfile
    with open('./inat_phen_MMRR_log.txt', 'w') as f:
        f.write(f'INAT MMRR LOG:\n{"-"*80}\n\n')

    # run analysis for each species with a non-unimodal hist
    # TODO: implement min number of samples here? (because actual n obs usually
    #       < expected n obs because of low positional_accuracy)
    mmrr_dict = {'tid': [],
                 'name': [],
                 'n': [],
                 'pct_miss': [],
                 'r2': [],
                 'int': [],
                 'int_t': [],
                 'int_p': [],
                 'geo': [],
                 'geo_t': [],
                 'geo_p': [],
                 'env': [],
                 'env_t': [],
                 'env_p': [],
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
                  'geo_dist': 'geo',
                  'geo_dist(t)': 'geo_t',
                  'geo_dist(p)': 'geo_p',
                  'env_dist': 'env',
                  'env_dist(t)': 'env_t',
                  'env_dist(p)': 'env_p',
                  'lsp_dist': 'lsp',
                  'lsp_dist(t)': 'lsp_t',
                  'lsp_dist(p)': 'lsp_p',
                  'F-statistic': 'f',
                  'F p-value': 'f_p',
                 }
    for i, row in res_gdf[res_gdf['signif_npeaks'] != 1].iterrows():
        start = time.time()
        try:
            tid = row['tid']
            name = row['name']
            print(f'\n\n{"-"*80}\n')
            print(f"running LSP MMRR on flowering observations for {name} (TID: {tid})...")
            # load observation points
            obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
            obs = gpd.read_file(os.path.join(obs_data_dir, obs_fn))
            res, n, pct_lsp_miss = fit_LSP_MMRR(obs,
                                                strict_coeffs=strict_coeffs,
                                                bioclim_nodata_val=bioclim_nodata_val,
                                                standardize_lsp_ts=True,
                                                interp_lsp_data=interp_lsp_data,
                                                neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                                                standardize_MMRR_data=True,
                                                fit_MMRR_intercept=True,
                                                MMRR_nperm=MMRR_nperm,
                                               )
            pprint.pprint(res)
            stop = time.time()
            runtime = stop-start
            runtime_msg = f"runtime: {np.round(runtime, 1)} s\n"
            print(runtime_msg)
            # log results
            pct_lsp_miss_msg = (f"{np.round(100*pct_lsp_miss, 2)}% "
                                "of sites are missing LSP coefficients")
            with open('./inat_phen_MMRR_log.txt', 'a') as f:
                f.write(f"{tid}: {name}\n")
                f.write(runtime_msg)
                f.write(f"{pct_lsp_miss_msg}\n\n")
                f.write(pprint.pformat(res))
                f.write(f'\n\n{"-"*80}\n\n')
        except Exception as e:
            res = {k: np.nan for k in label_dict.keys()}
            err_msg = f"\n\tERROR THROWN: {e}"
            print(f"{err_msg}\n\n\tmoving on...\n")
            # log results
            with open('./inat_phen_MMRR_log.txt', 'a') as f:
                f.write(f"{tid}: {name}\n")
                try:
                    f.write(f"{pct_lsp_miss_msg}\n\n")
                except Exception:
                    pass
                f.write(err_msg)
                f.write(f'\n\n{"-"*80}\n\n')
        mmrr_dict['tid'].append(tid)
        mmrr_dict['name'].append(name)
        # store results for this taxon
        try:
            mmrr_dict['n'].append(n)
        except Exception:
            mmrr_dict['n'].append(np.nan)
        try:
            mmrr_dict['pct_miss'].append(pct_lsp_miss)
        except Exception:
            mmrr_dict['pct_miss'].append(np.nan)
        for k in res.keys():
            mmrr_dict[label_dict[k]].append(res[k])

    # coerce results to pd.DataFrame
    mmrr_df = pd.DataFrame.from_dict(mmrr_dict)
    # sort by increasing P-values of interest and decreasing R^2 within that
    mmrr_df = mmrr_df.sort_values(by=['lsp_p', 'f_p', 'r2'],
                                  ascending=[True, True, False])
    mmrr_df.to_csv(mmrr_res_filename, index=False)

else:
    mmrr_df = pd.read_csv(mmrr_res_filename)


end_time = time.time()
runtime = end_time - start_time
print(f"\n\nTOTAL RUNTIME: {np.round(runtime/60/60, 2)} hours\n\n")

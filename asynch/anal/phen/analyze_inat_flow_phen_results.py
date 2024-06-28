#!/usr/bin/env python
# coding: utf-8

import os
import re
import h3
import sys
import pprint
import pyproj
import rasterstats
import numpy as np
import pandas as pd
import rasterio as rio
import geopandas as gpd
import rioxarray as rxr
import matplotlib as mpl
import cmcrameri.cm as cmc
import matplotlib.pyplot as plt
from shapely.geometry import box
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
# params and supporting data:
#------------------------------------------------------------------------------
# whether to make maps and whether to run MMRR models
make_peak_summary_maps = True
run_MMRRs = True

# whether to interpolate missing LSP coefficients, 
# and if so, how far away can rasterio look for neighbor to interpolate from?
# coefficient raster?
interp_sea_data = False
neigh_dist_sea_fill_tol = 2

# how many MMRR permutations to use
# NOTE: using this many permutations makes 0.001 the smallest possible
#       P-value, which is smaller than 0.05/number of taxa after filtering out
#       broadly distributed taxa, and thus allows for Bonferroni correction
MMRR_nperm = 999

# common equal-area projection to use
crs = 8857

# neighborhood size of asynchrony map to use
asynch_neigh_rad = 100

# whether to use the strict-masked (i.e., ag-removed) coefficients file for
# extracting fitted LSP patterns at flower observation localities
strict_coeffs = True

# minimum number of samples required to consider a taxa for the demonstrative
# plots
min_n_samps_for_demo_plots = 30

# directory where hex GeoJSON is to be saved
hex_data_dir = phf.EXTERNAL_INAT_DATA_DIR

# directory with downloaded iNat observation datasets
obs_data_dir = os.path.join(hex_data_dir, 'obs_data')

# load asynchrony map
asynch = rxr.open_rasterio(phf.ASYNCH_FILES[asynch_neigh_rad])[0]

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
    print("\n\treading hexes from file...\n")
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
    h3_gdf.to_file(os.path.join(hex_data_dir, peaks_hex_filename))
else:
    print("\n\treading hex-summarized phenology peak results from file...\n")
    h3_gdf = gpd.read_file(os.path.join(hex_data_dir, peaks_hex_filename))

if make_peak_summary_maps:
    fig = plt.figure(figsize=(9,16))
    gs = fig.add_gridspec(ncols=1, nrows=16)
    h3_gdf['prop_non1peak_signif'] = (h3_gdf['prop_0peak_signif'] +
                                      h3_gdf['prop_2pluspeak_signif'])
    cmap='viridis'
    # plot results for trees, non-trees, and all taxa combined
    res_cols = ['prop_non1peak_signif',
                'prop_0peak_signif',
                'prop_2pluspeak_signif',
               ]
    label_dict = {'prop_non1peak_signif': 'non-unimodal',
                  'prop_0peak_signif': 'no significant peaks',
                  'prop_2pluspeak_signif': '≥2 significant peaks',
                 }
    for i, res_col in enumerate(res_cols):
        gs_imin = 5*i
        gs_imax = 5*i+5+(1*(i==2))
        ax = fig.add_subplot(gs[gs_imin:gs_imax, 0])
        # add bottom axes for a colorbar
        if i == 2:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('bottom', size='7%', pad=0.2)
        # transform to equal-area projection and plot
        h3_gdf.to_crs(crs).plot(res_col,
                                cmap=cmap,
                                alpha=1,
                                zorder=0,
                                ax=ax,
                                edgecolor='white',
                                linewidth=0.2,
                                legend=False,
                                vmin=0,
                                #vmax=1,
                                vmax=np.nanpercentile(h3_gdf[res_col], 95),
                               )
        phf.plot_juris_bounds(ax,
                              lev1_alpha=0.6,
                              lev1_linewidth=0.05,
                              lev1_zorder=1,
                              lev0_linewidth=0.1,
                              lev0_zorder=2,
                              crs=crs,
                              strip_axes=True,
                              reset_axlims=False,
                             )
        ax.set_title(label_dict[res_col], fontdict={'fontsize':18})
        if i  == 2:
            scalcmap = plt.cm.ScalarMappable(cmap=cmap,
                                             norm=plt.Normalize(vmin=0, vmax=1),
                                            )
            plt.colorbar(scalcmap, cax=cax, orientation='horizontal')
            xticks = np.linspace(0, 1, 5)
            cax.set_xlabel('proportion of taxa',
                           fontdict={'fontsize': 15},
                          )
            cax.set_xticks(xticks, ['%0.2f' % t for t in xticks], size=9)
            cax.set_ylabel('')
            cax.set_yticks(())
    # adjust subplots and write to disk
    fig.subplots_adjust(top=0.96,
                        bottom=0.125,
                        right=0.96,
                        left=0.04,
                        hspace=0.9,
                       )
    fig.savefig(os.path.join(phf.FIGS_DIR,
                             'FIG_inat_nonunimodal_phen_summary_maps.png'),
                dpi=600,
               )


#------------------------------------------------------------------------------
# regress flowering observation doy distance on NIRv LSP data
#------------------------------------------------------------------------------
if not os.path.isfile('inat_flow_obs_LSP_MMRR_res.csv'):
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
        # get the lesser of the distance back to the earlier day of year or
        # forward to the same day of year next year
        d = sorted([doy1, doy2])
        dist = np.min((d[1]-d[0], d[0] + 365 - d[1]))
        return dist


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
            if strict_coeffs:
                coeffs_file = phf.COEFFS_STRICT_FILE
            else:
                coeffs_file = phf.COEFFS_FILE
            lsp_dist = phf.get_raster_info_points(coeffs_file,
                                                  pts,
                                                  'ts_pdm',
                                                  standardize=True,
                                                  fill_nans=interp_sea_data,
                                                  fill_tol=neigh_dist_sea_fill_tol,
                                                 )
            pct_miss = np.mean(pd.isnull(lsp_dist))
            pct_lsp_miss_msg = (f"{np.round(100*pct_miss, 2)} % "
                                "of lsp_dist consists of missing values")
            print(f"\n\t\t({pct_lsp_miss_msg})")
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
                       standardize=True,
                       intercept=True,
                       nperm=MMRR_nperm,
                      )
            pprint.pprint(res)
            # log results
            with open('./inat_phen_MMRR_log.txt', 'a') as f:
                f.write(f"{tid}: {name}\n")
                f.write(f"{pct_lsp_miss_msg}\n\n")
                f.write(pprint.pformat(res))
                f.write(f'\n\n{"-"*80}\n\n')
        except Exception as e:
            res = {k: np.nan for k in label_dict.keys()}
            err_msg = f"\n\tERROR THROWN: {e}\n\n\tmoving on...\n"
            print(f"{err_msg}\n\n\tmoving on...\n")
            # log results
            with open('./inat_phen_MMRR_log.txt', 'a') as f:
                f.write(f"{tid}: {name}\n")
                try:
                    f.write(f"{pct_lsp_miss_msg}\n\n")
                except Exception:
                    pass
                f.write(err_msg)
                f.write('\n\n{"-"*80}\n\n')
        mmrr_dict['tid'].append(tid)
        mmrr_dict['name'].append(name)
        # store results for this taxon
        try:
            mmrr_dict['n'].append(lsp_dist.shape[0])
        except Exception:
            mmrr_dict['n'].append(np.nan)
        try:
            mmrr_dict['pct_miss'].append(pct_miss)
        except Exception:
            mmrr_dict['pct_miss'].append(np.nan)
        for k in res.keys():
            mmrr_dict[label_dict[k]].append(res[k])

    # coerce results to pd.DataFrame
    mmrr_df = pd.DataFrame.from_dict(mmrr_dict)
    # sort by increasing P-values of interest and decreasing R^2 within that
    mmrr_df = mmrr_df.sort_values(by=['lsp_p', 'f_p', 'r2'],
                                  ascending=[True, True, False])
    mmrr_df.to_csv('inat_flow_obs_LSP_MMRR_res.csv', index=False)

else:
    mmrr_df = pd.read_csv('inat_flow_obs_LSP_MMRR_res.csv')


#------------------------------------------------------------------------------
# summarize MMRR results
#------------------------------------------------------------------------------

# what percent of non-unimodal flowering-histogram taxa have significant
# results for the flowering time~LSP correlation?
mmrr = mmrr_df[pd.notnull(mmrr_df['lsp_p'])]

# merge on alpha hulls and calculate hull bounding-box areas
mmrr_gdf = pd.merge(res_gdf.loc[:, ['tid', 'geometry']],
                    mmrr, on='tid',
                    how='right')
mmrr_gdf.loc[:, 'bbox_area_sqkm'] = [box(*g.bounds).area/(1000**2) for g in
                                     mmrr_gdf.to_crs(crs)['geometry']]
# move geom column to end again
mmrr_gdf = mmrr_gdf.iloc[:, [0]+[*range(2, mmrr_gdf.shape[1])]+[1]]

# subset to taxa with alpha hulls that are not evenly straddling the equator
# by at least 5 degrees on each side
# (i.e., simple way to filter out huge-range species that muddy the analysis)
spread_across_eq = []
for i, row in mmrr_gdf.iterrows():
    tid = row['tid']
    name = row['name']
    obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
    obs = gpd.read_file(os.path.join(obs_data_dir, obs_fn))
    min_y = np.min(obs.geometry.y)
    max_y = np.max(obs.geometry.y)
    spread_across = min_y <= 20 and max_y >= 20
    spread_across_eq.append(spread_across)
mmrr_gdf['spread_across_eq'] = spread_across_eq
mmrr_filt = mmrr_gdf[~mmrr_gdf['spread_across_eq']]
mmrr_filt = mmrr_filt.sort_values(['lsp_p'], ascending=[True])

# Bonferroni correction
alpha = 0.05/(len(mmrr_filt))
mmrr_filt.loc[:, 'lsp_p_sig_bonf_corr'] = mmrr_filt.loc[:, 'lsp_p']<=alpha

# print & log numbers & percentages of taxa passing filtering & analysis stages
n_nonunimodal_taxa = np.sum(res_gdf['signif_npeaks']!=1)
n_all_taxa = len(res_gdf)
sig_nonunimodal_msg = (f"{np.round(100*(n_nonunimodal_taxa/n_all_taxa))}% "
                       "of taxa are significantly non-unimodal "
                       f"({n_nonunimodal_taxa}/{n_all_taxa}).")
n_taxa_w_failed_MMRR = len(mmrr_df) - len(mmrr_gdf)
fail_fit_MMRR_msg = (f"{n_taxa_w_failed_MMRR} non-unimodal taxa failed "
                      "to fit a valid MMRR model.")
n_taxa_w_broad_distrs = len(mmrr_gdf) - len(mmrr_filt)
too_broad_distr_msg = (f"{n_taxa_w_broad_distrs} non-unimodal taxa dropped "
                        "because of broad, equator-crossing distributions.")
n_MMRR_taxa = len(mmrr_filt)
n_MMRR_taxa_LSP_signif = np.sum(mmrr_filt['lsp_p']<=0.05)
n_MMRR_taxa_LSP_signif_bonf_corr = np.sum(mmrr_filt['lsp_p_sig_bonf_corr'])
sig_MMRR_LSP_msg = (f"{np.round(100*(n_MMRR_taxa_LSP_signif/n_MMRR_taxa))}% "
                     "of non-unimodal taxa for which MMRRs were fitted "
                     "had significant LSP-distance coefficients (P<=0.05) "
                    f"({n_MMRR_taxa_LSP_signif}/{n_MMRR_taxa}).")
sig_MMRR_LSP_msg_bonf_corr = (
                f"{np.round(100*(n_MMRR_taxa_LSP_signif_bonf_corr/n_MMRR_taxa))}% "
                 "of non-unimodal taxa for which MMRRs were fitted "
                 "had significant LSP-distance coefficients after "
                 "Bonferroni correction "
                f"({n_MMRR_taxa_LSP_signif_bonf_corr}/{n_MMRR_taxa}).")
print(f"{fail_fit_hull_msg}\n\n")
print(f"{sig_nonunimodal_msg}\n\n")
print(f"{fail_fit_MMRR_msg}\n\n")
print(f"{too_broad_distr_msg}\n\n")
print(f"{sig_MMRR_LSP_msg}\n\n")
print(f"{sig_MMRR_LSP_msg_bonf_corr}\n\n")

with open('./inat_phen_MMRR_log.txt', 'a') as f:
    f.write(f'\n{"="*80}\n')
    f.write(f"{fail_fit_hull_msg}\n\n")
    f.write(f"{sig_nonunimodal_msg}\n\n")
    f.write(f"{fail_fit_MMRR_msg}\n\n")
    f.write(f"{too_broad_distr_msg}\n\n")
    f.write(f"{sig_MMRR_LSP_msg}\n\n")
    f.write(f"{sig_MMRR_LSP_msg_bonf_corr}\n\n")


# summarize MMRR results to hexes
mmrr_hex_filename = 'inat_hex_mmrr_results.json'
if not os.path.isfile(os.path.join(hex_data_dir, mmrr_hex_filename)):
    print("\n\tsummarizing MMRR results by hex...\n")
    h3_gdf = h3_gdf.loc[:, ['row_id', 'h3_id', 'x', 'y', 'geometry']]
    # summarize results within hexes
    new_cols = {'n_taxa': [],
                'mean_n': [],
                'mean_pct_miss': [],
                'mean_lsp_p': [],
                'mean_lsp': [],
                'prop_lsp_signif': [],
                'prop_lsp_signif_0p05': [],
                'mean_geo_p': [],
                'mean_geo': [],
                'prop_geo_signif': [],
                'mean_int_p': [],
                'mean_int': [],
                'prop_int_signif': [],
                'mean_f_p': [],
                'prop_mod_signif': [],
                'mean_r2': [],
               }
    for i, row in h3_gdf.iterrows():
        # get all taxa that intersect with this hex
        gdf_intsxn = mmrr_filt[mmrr_filt.intersects(row['geometry'])]
        new_cols['n_taxa'].append(len(gdf_intsxn))
        new_cols['mean_n'].append(np.mean(gdf_intsxn['n']))
        new_cols['mean_pct_miss'].append(np.mean(gdf_intsxn['pct_miss']))
        new_cols['mean_lsp_p'].append(np.mean(gdf_intsxn['lsp_p']))
        new_cols['mean_lsp'].append(np.mean(gdf_intsxn['lsp']))
        new_cols['prop_lsp_signif'].append(np.mean(gdf_intsxn['lsp_p']<=alpha))
        new_cols['prop_lsp_signif_0p05'].append(np.mean(gdf_intsxn['lsp_p']<=0.05))
        new_cols['mean_geo_p'].append(np.mean(gdf_intsxn['geo_p']))
        new_cols['mean_geo'].append(np.mean(gdf_intsxn['geo']))
        new_cols['prop_geo_signif'].append(np.mean(gdf_intsxn['geo_p']<=alpha))
        new_cols['mean_int_p'].append(np.mean(gdf_intsxn['int_p']))
        new_cols['mean_int'].append(np.mean(gdf_intsxn['int']))
        new_cols['prop_int_signif'].append(np.mean(gdf_intsxn['int_p']<=alpha))
        new_cols['mean_f_p'].append(np.mean(gdf_intsxn['f_p']))
        new_cols['prop_mod_signif'].append(np.mean(gdf_intsxn['f_p']<=alpha))
        new_cols['mean_r2'].append(np.mean(gdf_intsxn['r2']))
    for col, vals in new_cols.items():
        h3_gdf[col] = vals
        if col.startswith('prop'):
            assert (np.all(h3_gdf[pd.notnull(h3_gdf[col])][col]>=0) and
                    np.all(h3_gdf[pd.notnull(h3_gdf[col])][col]<=1))
    h3_gdf.to_file(os.path.join(hex_data_dir, mmrr_hex_filename))
else:
    print("\n\treading hex-summarized LSP MMRR results from file...\n")
    h3_gdf = gpd.read_file(os.path.join(hex_data_dir, mmrr_hex_filename))
# map hex-summarized results
fig = plt.figure(figsize=(9,5))
ax = fig.add_axes((0.01, 0.01, 0.97, 0.97))
cmap='viridis'
# add bottom axes for a colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes('bottom', size='7%', pad=0.2)
# transform to equal-area projection and plot
h3_gdf.to_crs(crs).plot('prop_lsp_signif_0p05',
                        cmap=cmap,
                        alpha=1,
                        zorder=2,
                        ax=ax,
                        edgecolor='white',
                        linewidth=0.2,
                        legend=False,
                        vmin=0,
                        vmax=1,
                        #vmax=np.nanpercentile(h3_gdf['prop_lsp_signif'], 90),
                       )
phf.plot_juris_bounds(ax,
                      lev1_alpha=0.6,
                      lev1_linewidth=0.05,
                      lev1_zorder=0,
                      lev0_linewidth=0.1,
                      lev0_zorder=1,
                      crs=crs,
                      strip_axes=True,
                      reset_axlims=False,
                     )
scalcmap = plt.cm.ScalarMappable(cmap=cmap,
                                 norm=plt.Normalize(vmin=0, vmax=1),
                                )
plt.colorbar(scalcmap, cax=cax, orientation='horizontal')
xticks = np.linspace(0, 1, 5)
cax.set_xlabel('proportion taxa with significant LSP MMRR coefficient',
               fontdict={'fontsize': 11},
              )
cax.set_xticks(xticks, ['%0.2f' % t for t in xticks], size=9)
cax.set_ylabel('')
cax.set_yticks(())
fig.show()


#------------------------------------------------------------------------------
# plot demonstrative taxa
#------------------------------------------------------------------------------

# drop small sample sizes (i.e., where too many LSP pixels were masked out)
mmrr_filt = mmrr_filt[mmrr_filt['n'] >=min_n_samps_for_demo_plots]

# get just the top LSP-significant taxa
mmrr_filt = mmrr_filt.sort_values(by=['lsp_p', 'f_p', 'r2'],
                                  ascending=[True, True, False],
                                 )

# find example taxa in areas of interest:
countries = gpd.read_file(phf.ADM0_BOUNDS).to_crs(crs).loc[:, ['adm0_a3', 'geometry']]
subnational = gpd.read_file(phf.SELECT_ADM1_BOUNDS).to_crs(crs)

# Baja
baja_n = subnational[subnational['NAME_1'] == 'BajaCalifornia']
baja_s = subnational[subnational['NAME_1'] == 'BajaCaliforniaSur']
baja_n_taxa = baja_n.loc[:, ['geometry',
                             'COUNTRY']].sjoin(mmrr_filt.to_crs(crs))
baja_s_taxa = baja_s.loc[:, ['geometry',
                             'COUNTRY']].sjoin(mmrr_filt.to_crs(crs))
baja_taxa = baja_n_taxa[baja_n_taxa['tid'].isin(baja_s_taxa['tid'].values)]

# Texas
tx = subnational[subnational['NAME_1'] == 'Texas']
tx_taxa = tx.loc[:, ['geometry',
                     'COUNTRY']].sjoin(mmrr_filt.to_crs(crs))

# northen Andes
andes= countries[countries['adm0_a3'].isin(['COL',
                                               'ECU',])].dissolve()
andes_taxa = andes.loc[:, ['geometry', 'adm0_a3']].sjoin(
                                                    mmrr_filt.to_crs(crs))

# Eastern Brazil
n_e_brz = subnational[subnational['NAME_1'].isin(['Bahia', 'Sergipe',
                                                 'Alagoas', 'Pernambuco',
                                                 'RioGrandedoNorte',
                                                 'Paraíba'])].dissolve()
s_e_brz = subnational[subnational['NAME_1'].isin(['EspíritoSanto',
                                                 'RiodeJaneiro',
                                                 'MinasGerais'])].dissolve()
n_e_brz_taxa = n_e_brz.loc[:, ['geometry', 'COUNTRY']].sjoin(mmrr_filt.to_crs(crs))
s_e_brz_taxa = s_e_brz.loc[:, ['geometry', 'COUNTRY']].sjoin(mmrr_filt.to_crs(crs))
e_brz_taxa = n_e_brz_taxa[n_e_brz_taxa['tid'].isin(s_e_brz_taxa['tid'].values)]

# South Africa
zaf = countries[countries['adm0_a3'] == 'ZAF'].dissolve()
zaf_taxa = zaf.loc[:, ['geometry', 'adm0_a3']].sjoin(
                                                    mmrr_filt.to_crs(crs))

# Australia
aus = countries[countries['adm0_a3'] == 'AUS'].dissolve()
aus_taxa = aus.loc[:, ['geometry', 'adm0_a3']].sjoin(
                                                    mmrr_filt.to_crs(crs))

# re-sort main and focal dfs
for df in [baja_taxa, tx_taxa, andes_taxa, e_brz_taxa, zaf_taxa, aus_taxa]:
    df.sort_values('lsp_p', ascending=True, inplace=True)


# load the scaled, ITCZ-folded EOFs, for data viz of phen significant MMRR results
dataset = 'NIRv'
normts_file_substr = '_normts'
mask_filename_ext = ''
eofs = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
            '%s_4_EOFs_sqrt_coswts%s_SCALED_FOLDED_EPSG-8857.tif') %(dataset,
                                                      normts_file_substr))[:3]
if strict_coeffs:
    strict_coeffs_xarr = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                                     phf.COEFFS_STRICT_FILE),
                                        masked=True)[0].rio.reproject(crs)
    eofs = eofs.where(pd.notnull(strict_coeffs_xarr), np.nan)



def plot_annual_flowering_clock(obs, eofs, ax):
    for i, row in obs.iterrows():
        g = row['geometry']
        e = eofs.sel(x=g.x, y=g.y, method='nearest').values
        if np.all(pd.notnull(e)):
            c = tuple(e)
            d = row['doy_circ']
            # add jitter
            d = d + np.random.normal(0, 0.1)
            x = np.sin(d)
            y = np.cos(d)
            xs = [0, x]
            ys = [0, y]
            #ax.scatter(x, y, c=c, s=50, edgecolor='black', alpha=0.9)
            ax.plot(xs, ys, color=c, linewidth=1, alpha=0.9)
    circle_vals = np.arange(0, 2*np.pi+0.01, 0.001)
    circle_xs = np.sin(circle_vals)
    circle_ys = np.cos(circle_vals)
    months = ['Jan', 'Apr', 'Jul', 'Oct']
    # NOTE: slightly offset so as to center text
    month_xs = [-0.05, 1.02, -0.05, -1.15]
    month_ys = [1.05, 0, -1.1, 0]
    for m, x, y in zip(months, month_xs, month_ys):
        ax.text(x, y, m, size=10)
    ax.plot(circle_xs, circle_ys, '-k', linewidth=0.5)
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.set_aspect('equal')


def determine_map_bounds(obs):
    coords = [np.min(obs.geometry.x),
              np.min(obs.geometry.y),
              np.max(obs.geometry.x),
              np.max(obs.geometry.y),
             ]
    ranges = [np.abs(coords[2] - coords[0]),
              np.abs(coords[3] - coords[1]),
             ]
    margins = [-0.05*ranges[0],
               -0.05*ranges[1],
                0.05*ranges[0],
                0.05*ranges[1],
              ]
    # NOTE: assuming no coords are exactly 0
    bounds = [coords[i] + margins[i] for i in range(4)]
    return bounds


def plot_signif_taxon(tid,
                      mmrr_df,
                      scaled_eofs,
                      strict_coeffs=strict_coeffs,
                      fill_nans=False,
                      fill_tol=None,
                      mmrr_res=None,
                      map=True,
                      savefig=False,
                     ):
    if fill_nans:
        assert fill_tol is not None and isinstance(fill_tol, int)
    if mmrr_res is not None:
        assert ((isinstance(mmrr_res, dict) or
                 isinstance(mmrr_res, OrderedDict)) and
                'lsp_dist' in mmrr_res.keys() and
                'lsp_dist(p)' in mmrr_res.keys())
    name = mmrr_df[mmrr_df['tid'] == tid]['name'].values[0]
    pval = mmrr_df[mmrr_df['tid'] == tid]['lsp_p'].values[0]
    obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
    obs = gpd.read_file(os.path.join(obs_data_dir, obs_fn))
    pts = np.vstack([obs.geometry.x, obs.geometry.y]).T
    obs = obs.to_crs(crs)
    if strict_coeffs:
        coeffs_file = phf.COEFFS_STRICT_FILE
    else:
        coeffs_file = phf.COEFFS_FILE
    lsp_ts = phf.get_raster_info_points(coeffs_file,
                                        pts,
                                        'ts',
                                        standardize=True,
                                        fill_nans=fill_nans,
                                        fill_tol=fill_tol,
                                       )
    missing = np.any(pd.isnull(lsp_ts), axis=1)
    print((f"\n{np.round(100*(np.sum(missing)/lsp_ts.shape[0]))}% of sites "
            "have no LSP data\n"))
    fig = plt.figure()
    ax = fig.add_axes([0.04, 0.54, 0.9, 0.4])
    doy_scale = (obs['doy'].values - np.min(obs['doy']))/(np.max(obs['doy']) -
                                                          np.min(obs['doy']))
    doy_colors = cmc.brocO_r(doy_scale)
    for i in range(len(obs)):
        if not missing[i]:
            color = doy_colors[i]
            ax.plot(lsp_ts[i, :],
                    color=color,
                    alpha=0.7,
                    linewidth=1.6,
                    zorder=0,
                   )
            dot_color = np.array(scaled_eofs.sel(x=obs.iloc[i, :]['geometry'].x,
                                              y=obs.iloc[i, :]['geometry'].y,
                                              method='nearest',
                                             )).reshape((1, -1))
            d = obs.iloc[i, :]['doy']
            assert d <=366
            ax.scatter(x=np.clip(d, a_min=None, a_max=365),
                       y=lsp_ts[i, np.clip(int(d), a_min=None, a_max=364)],
                       c=dot_color,
                       s=50,
                       edgecolor='black',
                       linewidth=0.4,
                       zorder=1,
                      )
    ax.set_xlim(0, 364)
    title = f"{name} (TID: {tid})"
    title = title + " ($P_{LSP}=%0.2f$)" % pval
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
    cbar = ColorbarBase(cax, cmap=cmc.brocO_r, norm=norm, orientation='vertical')
    cbar.set_label('flowering observation day of year', fontsize=12)
    ax = fig.add_axes([0.04, 0.04, 0.4, 0.4])
    ax.set_aspect('equal')
    plot_annual_flowering_clock(obs, eofs, ax)
    if map:
        ax = fig.add_axes([0.5, 0.04, 0.42, 0.42])
        phf.plot_juris_bounds(ax,
                              lev1=False,
                              lev0_linewidth=0.1,
                              lev0_zorder=2,
                              crs=obs.crs,
                              strip_axes=True,
                             )
        obs.plot('doy',
                 cmap=cmc.brocO_r,
                 alpha=0.8,
                 zorder=3,
                 ax=ax,
                 vmin=0,
                 vmax=365,
                 edgecolor='black',
                 linewidth=0.4,
                )
        eofs.plot.imshow(ax=ax,
                         zorder=0,
                        )
        bounds = determine_map_bounds(obs)
        ax.set_xlim(bounds[0], bounds[2])
        ax.set_ylim(bounds[1], bounds[3])
        ax.set_aspect('equal')
        phf.strip_axes_labels_and_ticks(ax)
    if savefig:
        fig.savefig(os.path.join(phf.FIGS_DIR,
                                 f"MMRR_res_figs/TID_{tid}_{name.replace(' ', '_')}_top_MMRR_res.png"),
                    dpi=600,
                   )
    return fig




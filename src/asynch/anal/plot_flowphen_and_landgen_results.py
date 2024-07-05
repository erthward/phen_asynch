#!/usr/bin/env python
# coding: utf-8














"""
 TODO:
 - finalize taxa
 - refine map bounds
 - coordinate cluster colors for multiple taxa within a region
 - move colorbar to top of global map and improve label
 - plot scree plot for each taxon's lsp clusters and double-check k! (see marsypianthes)
 - plot demo taxa for sw and saf
 - arrange and clean up draft plot
 - figure out month labels for radar plots, and brainstorm improvements for clarity...
 - change colored scatterplots to pie charts?
 - could clustering be used to determine folding points in eof map???
"""










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
from shapely.geometry import Point
from collections import OrderedDict
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from shapely.geometry import Polygon, MultiPolygon
from mpl_toolkits.axes_grid1 import make_axes_locatable

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf
from MMRR import MMRR

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/asynch/anal/gen/xiphorhynchus/'))
import test_xiphorhynchus_fuscus as xf

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/asynch/anal/gen/rhinella/'))
import test_rhinella_granulosa as rg



#------------------------------------------------------------------------------
# set overarching params and load supporting data:
#------------------------------------------------------------------------------

# whether to interpolate missing LSP coefficients, 
# and if so, how far away can rasterio look for neighbor to interpolate from?
# coefficient raster?
interp_lsp_data = False
neigh_dist_lsp_fill_tol = 2

# common equal-area projection to use
crs = 8857

# whether to use the strict-masked (i.e., ag-removed) coefficients file for
# extracting fitted LSP patterns at flower observation localities
strict_coeffs = True
if strict_coeffs:
    coeffs_file = phf.COEFFS_STRICT_FILE
else:
    coeffs_file = phf.COEFFS_FILE

# minimum number of iNat flowering observation samples required to consider a taxon for plotting
min_n_inat_samps_for_demo_plots = 30

# directory where iNat hex GeoJSONs are to be saved
inat_hex_data_dir = phf.EXTERNAL_INAT_DATA_DIR

# directory with downloaded iNat observation datasets
inat_obs_data_dir = os.path.join(inat_hex_data_dir, 'obs_data')

# load the scaled, ITCZ-folded EOFs, for data viz of phen significant MMRR results
eofs = rxr.open_rasterio(phf.EOFS_PREPPED_FILE)[:3]
if strict_coeffs:
    strict_coeffs_xarr = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                                        coeffs_file),
                                           masked=True)[0].rio.reproject(crs)
    eofs = eofs.where(pd.notnull(strict_coeffs_xarr), np.nan)


#------------------------------------------------------------------------------
# summarize and plot flower-phenology analysis results
#------------------------------------------------------------------------------

# create the overall figure
fig = plt.figure(figsize=(16,22))
gs = fig.add_gridspec(ncols=180, nrows=220)

# create the iNat peak-analysis supp figure
fig_supp = plt.figure(figsize=(9,16))

# plot inat flower phenology peak-analysis results
peaks_hex_filename = 'inat_hex_results.json'
inat_h3_gdf = gpd.read_file(os.path.join(inat_hex_data_dir, peaks_hex_filename))

inat_h3_gdf['prop_non1peak_signif'] = (inat_h3_gdf['prop_0peak_signif'] +
                                  inat_h3_gdf['prop_2pluspeak_signif'])
cmap='viridis'
# plot results for trees, non-trees, and all taxa combined
res_cols = ['prop_non1peak_signif',
            'prop_0peak_signif',
            'prop_2pluspeak_signif',
            'prop_non1peak_signif',
           ]
label_dict = {'prop_non1peak_signif': 'non-unimodal',
              'prop_0peak_signif': 'no significant peaks',
              'prop_2pluspeak_signif': '≥2 significant peaks',
             }
for i, res_col in enumerate(res_cols):
    gs_imin = 5*i
    gs_imax = 5*i+5+(1*(i==2))
    if i == 3:
        ax = fig.add_subplot(gs[:45, :70])
    else:
        ax = fig_supp.add_subplot(gs[gs_imin:gs_imax, 0])
    # add bottom axes for a colorbar
    if i >= 2:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('bottom', size='7%', pad=0.2)
    # transform to equal-area projection and plot
    inat_h3_gdf.to_crs(crs).plot(res_col,
                            cmap=cmap,
                            alpha=1,
                            zorder=0,
                            ax=ax,
                            edgecolor='white',
                            linewidth=0.2,
                            legend=False,
                            vmin=0,
                            #vmax=1,
                            vmax=np.nanpercentile(inat_h3_gdf[res_col], 95),
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
    if i < 3:
        ax.set_title(label_dict[res_col], fontdict={'fontsize':18})
    if i  >= 2:
        scalcmap = plt.cm.ScalarMappable(cmap=cmap,
                                         norm=plt.Normalize(vmin=0, vmax=1),
                                        )
        plt.colorbar(scalcmap, cax=cax, orientation='horizontal')
        xticks = np.linspace(0, 1, 5)
        cax.set_xlabel('proportion of taxa present',
                       fontdict={'fontsize': 15},
                      )
        cax.set_xticks(xticks, ['%0.2f' % t for t in xticks], size=9)
        cax.set_ylabel('')
        cax.set_yticks(())

# save the supplemental figure
fig_supp.subplots_adjust(top=0.96,
                         bottom=0.125,
                         right=0.96,
                         left=0.04,
                         hspace=0.9,
                        )
fig_supp.savefig(os.path.join(phf.FIGS_DIR,
                              'FIG_SUPP_inat_taxa_npeaks_summary_maps.png'),
                 dpi=600,
                )


#------------------------------------------------------------------------------
# summarize and plot flowering_phenology MMRR results
#------------------------------------------------------------------------------
mmrr_res_filename = './phen/inat_flower_doy_LSP_MMRR_res.csv'
inat_mmrr_df = pd.read_csv(mmrr_res_filename)

# drop all taxa for which MMRR models could not be fit
# (typically, because too few samples coincided with non-masked LSP pixels)
inat_mmrr_df_valid_pval = inat_mmrr_df[pd.notnull(inat_mmrr_df['lsp_p'])]

# merge on alpha hulls and calculate hull bounding-box areas
inat_res_gdf = gpd.read_file('./phen/inat_flower_phen_results.json')
len_b4 = len(inat_res_gdf)
inat_res_gdf = inat_res_gdf[pd.notnull(inat_res_gdf['geometry'])]
len_af = len(inat_res_gdf)
fail_fit_hull_msg = (f"\n\n{len_b4-len_af} taxa dropped because of "
                     "failure to fit observation-range alpha hulls.")
inat_mmrr_gdf = pd.merge(inat_res_gdf.loc[:, ['tid', 'geometry']],
                    inat_mmrr_df_valid_pval, on='tid',
                    how='right')
# calculate a column of the fitted number of histogram peaks,
# based on the KDE npeaks and its P-value
sig_npeaks = [row['npeaks'] if row['npeaks_pval']<=(0.05/len(inat_res_gdf)) else
                            0 for i, row in inat_res_gdf.iterrows()]
inat_res_gdf.loc[:, 'signif_npeaks'] = sig_npeaks

# move geom column to end again
inat_mmrr_gdf = inat_mmrr_gdf.iloc[:, [0]+[*range(2, inat_mmrr_gdf.shape[1])]+[1]]

# flag and filter out taxa with alpha hulls that cross
# the equator by 10 degrees or more on each side
# (a simple way to filter out huge-range species that muddy the analysis)
extreme_lat_range = []
for i, row in inat_mmrr_gdf.iterrows():
    tid = row['tid']
    name = row['name']
    obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
    obs = gpd.read_file(os.path.join(inat_obs_data_dir, obs_fn))
    min_y = np.min(obs.geometry.y)
    max_y = np.max(obs.geometry.y)
    has_extreme_lat_range = min_y <= -10 and max_y >= 10
    extreme_lat_range.append(has_extreme_lat_range)
inat_mmrr_gdf['extreme_lat_range'] = extreme_lat_range

# Bonferroni correction
# (using only the taxa that don't have extreme latitudinal ranges)
p_bonf_corr = 0.05/(np.sum(np.invert(inat_mmrr_gdf['extreme_lat_range'])))
inat_mmrr_gdf.loc[:, 'lsp_p_sig_bonf_corr'] = inat_mmrr_gdf.loc[:, 'lsp_p']<=p_bonf_corr
inat_mmrr_gdf.loc[inat_mmrr_gdf['extreme_lat_range'], 'lsp_p_sig_bonf_corr'] = np.nan

# sort rows by significance and then by sample size
inat_mmrr_gdf = inat_mmrr_gdf.sort_values(by=['extreme_lat_range', 'lsp_p', 'n'],
                                          ascending=[True, True, False],
                                         )

# save MMRR results formatted as a supplemental table
inat_mmrr_gdf_for_table = inat_mmrr_gdf.loc[:, ['name',
                                                'tid',
                                                'pct_miss',
                                                'n',
                                                'int',
                                                'int_t',
                                                'int_p',
                                                'geo',
                                                'geo_t',
                                                'geo_p',
                                                'env',
                                                'env_t',
                                                'env_p',
                                                'lsp',
                                                'lsp_t',
                                                'lsp_p',
                                                'lsp_p_sig_bonf_corr',
                                                'f',
                                                'f_p',
                                                'r2',
                                               ]]
inat_mmrr_gdf_for_table.columns = ['taxon name',
                                   'taxon iNat ID',
                                   '% observation locations missing LSP data',
                                   'MMRR n',
                                   'MMRR intercept',
                                   'MMRR intercept t-stat',
                                   'MMRR intercept P-value',
                                   'MMRR geographic distance coefficient',
                                   'MMRR geographic distance coefficient t-stat',
                                   'MMRR geographic distance coefficient P-value',
                                   'MMRR environmental distance coefficient',
                                   'MMRR environmental distance coefficient t-stat',
                                   'MMRR environmental distance coefficient P-value',
                                   'MMRR LSP distance coefficient',
                                   'MMRR LSP distance coefficient t-stat',
                                   'MMRR LSP distance coefficient P-value',
                                   'MMRR LSP distance coefficient signif. after Bonferroni correction?',
                                   'MMRR F-stat',
                                   'MMRR F-stat P-value',
                                   'MMRR R2',
                                  ]
inat_mmrr_gdf_for_table.to_csv(os.path.join(phf.TABS_DIR,
                                            'TAB_SUPP_iNat_MMRR_results.csv'),
                               index=False,
                              )

# drop extreme-latitudinal-range taxa from further analysis
inat_mmrr_filt = inat_mmrr_gdf[~inat_mmrr_gdf['extreme_lat_range']]

# print & log numbers & percentages of taxa passing the different filtering & analysis stages
n_nonunimodal_taxa = np.sum(inat_res_gdf['signif_npeaks']!=1)
n_all_taxa = len(inat_res_gdf)
sig_nonunimodal_msg = (f"{np.round(100*(n_nonunimodal_taxa/n_all_taxa))}% "
                       "of taxa are significantly non-unimodal "
                       f"({n_nonunimodal_taxa}/{n_all_taxa}).")
n_taxa_w_failed_MMRR = len(inat_mmrr_df) - len(inat_mmrr_df_valid_pval)
fail_fit_MMRR_msg = (f"{n_taxa_w_failed_MMRR} non-unimodal taxa failed "
                      "to fit a valid MMRR model.")
n_taxa_w_broad_distrs = len(inat_mmrr_gdf) - len(inat_mmrr_filt)
too_broad_distr_msg = (f"{n_taxa_w_broad_distrs} non-unimodal taxa dropped "
                        "because of broad, equator-crossing distributions.")
n_MMRR_taxa = len(inat_mmrr_filt)
n_MMRR_taxa_LSP_signif_0p05 = np.sum(inat_mmrr_filt['lsp_p']<=0.05)
n_MMRR_taxa_LSP_signif_bonf_corr = np.sum(inat_mmrr_filt['lsp_p_sig_bonf_corr'])
sig_MMRR_LSP_msg = (f"{np.round(100*(n_MMRR_taxa_LSP_signif_0p05/n_MMRR_taxa))}% "
                     "of non-unimodal taxa for which MMRRs were fitted "
                     "had significant LSP-distance coefficients (P<=0.05) "
                    f"({n_MMRR_taxa_LSP_signif_0p05}/{n_MMRR_taxa}).")
sig_MMRR_LSP_msg_bonf_corr = (
                f"{np.round(100*(n_MMRR_taxa_LSP_signif_bonf_corr/n_MMRR_taxa))}% "
                 "of non-unimodal taxa for which MMRRs were fitted "
                 "had significant LSP-distance coefficients after "
                f"Bonferroni correction (P-value <={np.round(p_bonf_corr, 5)}) "
                f"({n_MMRR_taxa_LSP_signif_bonf_corr}/{n_MMRR_taxa}).")
print(f"{fail_fit_hull_msg}\n\n")
print(f"{sig_nonunimodal_msg}\n\n")
print(f"{fail_fit_MMRR_msg}\n\n")
print(f"{too_broad_distr_msg}\n\n")
print(f"{sig_MMRR_LSP_msg}\n\n")
print(f"{sig_MMRR_LSP_msg_bonf_corr}\n\n")

# log the summary info, if not already logged
MMRR_sum_logged = False
with open('./phen/inat_phen_MMRR_log.txt', 'r') as f:
    if re.search('={80}', f.read()) is not None:
        MMRR_sum_logged = True
if not MMRR_sum_logged:
    with open('./phen/inat_phen_MMRR_log.txt', 'a') as f:
        f.write(f'\n{"="*80}\n')
        f.write(f"{fail_fit_hull_msg}\n\n")
        f.write(f"{sig_nonunimodal_msg}\n\n")
        f.write(f"{fail_fit_MMRR_msg}\n\n")
        f.write(f"{too_broad_distr_msg}\n\n")
        f.write(f"{sig_MMRR_LSP_msg}\n\n")
        f.write(f"{sig_MMRR_LSP_msg_bonf_corr}\n\n")

# summarize MMRR results to hexes
mmrr_hex_filename = 'inat_hex_MMRR_results.json'
if not os.path.isfile(os.path.join(inat_hex_data_dir, mmrr_hex_filename)):
    print("\n\tsummarizing MMRR results by hex...\n")
    inat_h3_gdf = inat_h3_gdf.loc[:, ['row_id', 'h3_id', 'x', 'y', 'geometry']]
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
    for i, row in inat_h3_gdf.iterrows():
        # get all taxa that intersect with this hex
        gdf_intsxn = inat_mmrr_filt[inat_mmrr_filt.intersects(row['geometry'])]
        new_cols['n_taxa'].append(len(gdf_intsxn))
        new_cols['mean_n'].append(np.mean(gdf_intsxn['n']))
        new_cols['mean_pct_miss'].append(np.mean(gdf_intsxn['pct_miss']))
        new_cols['mean_lsp_p'].append(np.mean(gdf_intsxn['lsp_p']))
        new_cols['mean_lsp'].append(np.mean(gdf_intsxn['lsp']))
        new_cols['prop_lsp_signif'].append(np.mean(gdf_intsxn['lsp_p']<=p_bonf_corr))
        new_cols['prop_lsp_signif_0p05'].append(np.mean(gdf_intsxn['lsp_p']<=0.05))
        new_cols['mean_geo_p'].append(np.mean(gdf_intsxn['geo_p']))
        new_cols['mean_geo'].append(np.mean(gdf_intsxn['geo']))
        new_cols['prop_geo_signif'].append(np.mean(gdf_intsxn['geo_p']<=p_bonf_corr))
        new_cols['mean_int_p'].append(np.mean(gdf_intsxn['int_p']))
        new_cols['mean_int'].append(np.mean(gdf_intsxn['int']))
        new_cols['prop_int_signif'].append(np.mean(gdf_intsxn['int_p']<=p_bonf_corr))
        new_cols['mean_f_p'].append(np.mean(gdf_intsxn['f_p']))
        new_cols['prop_mod_signif'].append(np.mean(gdf_intsxn['f_p']<=p_bonf_corr))
        new_cols['mean_r2'].append(np.mean(gdf_intsxn['r2']))
    for col, vals in new_cols.items():
        inat_h3_gdf[col] = vals
        if col.startswith('prop'):
            assert (np.all(inat_h3_gdf[pd.notnull(inat_h3_gdf[col])][col]>=0) and
                    np.all(inat_h3_gdf[pd.notnull(inat_h3_gdf[col])][col]<=1))
    inat_h3_gdf.to_file(os.path.join(inat_hex_data_dir, mmrr_hex_filename))
else:
    print("\n\treading hex-summarized LSP MMRR results from file...\n")
    inat_h3_gdf = gpd.read_file(os.path.join(inat_hex_data_dir, mmrr_hex_filename))


#------------------------------------------------------------------------------
# plot demonstrative taxa
#------------------------------------------------------------------------------

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


def plot_sig_taxon(inat_mmrr_res_df_row,
                   colors=np.array(['#2d5098', '#ca1957']),
                   subset_poly=None,
                   subset_poly_buffer=None,
                   save_it=False,
                  ):
    # read observations
    tid = inat_mmrr_res_df_row['tid']
    name = inat_mmrr_res_df_row['name']
    obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
    obs = gpd.read_file(os.path.join(inat_obs_data_dir, obs_fn))
    if subset_poly is not None:
        if subset_poly_buffer is not None:
            poly = subset_poly.buffer(subset_poly_buffer)
        else:
            poly = subset_poly
        obs = obs[obs.to_crs(poly.crs).within(poly.iloc[0])]
    bounds = determine_map_bounds(obs)
    map_xlim = bounds[0::2]
    map_ylim = bounds[1::2]
    map_lim_df = pd.DataFrame({'geometry':[Point(*lims) for lims in zip(map_xlim,
                                                                    map_ylim)],
                               'idx': range(2)})
    map_lim_gdf = gpd.GeoDataFrame(map_lim_df, geometry='geometry', crs=4326)
    map_lim_proj_coords = map_lim_gdf.to_crs(crs).get_coordinates().values
    map_xlim = map_lim_proj_coords[:, 0]
    map_ylim = map_lim_proj_coords[:, 1]
    fig = plt.figure(figsize=(10,5))
    gs = fig.add_gridspec(ncols=100, nrows=50)
    ax_rgb = fig.add_subplot(gs[:, :50])
    eofs.plot.imshow(ax=ax_rgb)
    obs.to_crs(crs).plot(color='black', markersize=15, alpha=0.7, ax=ax_rgb)
    phf.plot_juris_bounds(ax_rgb,
                          crs=crs,
                          strip_axes=True,
                          reset_axlims=False,
                         )
    ax_rgb.set_xlim(map_xlim)
    ax_rgb.set_ylim(map_ylim)
    ax_map = fig.add_subplot(gs[:15, 50:])
    ax_ts = fig.add_subplot(gs[15:30, 50:])
    ax_radar = fig.add_subplot(gs[30:50, 50:], projection='polar')
    phf.plot_flowerdate_LSP_comparison(flower_obs=obs,
                                       ax_map=ax_map,
                                       ax_ts=ax_ts,
                                       ax_radar=ax_radar,
                                       plot_crs=crs,
                                       map_xlim=map_xlim,
                                       map_ylim=map_ylim,
                                       interp_lsp_data=interp_lsp_data,
                                       neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                                       radar_alpha=0.5,
                                       radar_width_shrink_factor=0.9,
                                      )
    if save_it:
        fig.savefig(os.path.join(phf.FIGS_DIR,
                             'MMRR_res_figs',
            f"TID_{tid}_{name.replace(' ', '_')}_MMRR_res_map_ts_radar.png"))
    else:
        fig.show()
        input('<Enter> to close...')
    plt.close('all')


# drop small sample sizes (i.e., where too many LSP pixels were masked out)
inat_mmrr_filt_adeq_n = inat_mmrr_filt[inat_mmrr_filt['n'] >=min_n_inat_samps_for_demo_plots]

# set cluster colors
clust_colors=np.array(['#2d5098', # blue
                       '#ca1957', # red
                       '#baeb34', # lime
                       '#8b54bf', # purple
                      ])

# select and plot exemplary taxa in North American southwest
sw_names = ['Datura wrightii',
            'Allionia incarnata',
            'Xanthisma spinulosum',
           ]
sw_taxa = inat_mmrr_filt_adeq_n[inat_mmrr_filt_adeq_n['name'].isin(sw_names)]
# set values of K (based on manually inspected scree plots of LSP clustering)
K_vals = {'Datura wrightii': 3,
          'Allionia incarnata': 3,
          'Xanthisma spinulosum': 2,
         }
# set map bounding box
map_xlim = [-126, -94]
map_ylim = [21, 45]
map_lim_df = pd.DataFrame({'geometry':[Point(*lims) for lims in zip(map_xlim,
                                                                    map_ylim)],
                           'idx': range(2)})
map_lim_gdf = gpd.GeoDataFrame(map_lim_df, geometry='geometry', crs=4326)
map_lim_proj_coords = map_lim_gdf.to_crs(crs).get_coordinates().values
map_xlim = map_lim_proj_coords[:, 0]
map_ylim = map_lim_proj_coords[:, 1]
ax_rgb_map = fig.add_subplot(gs[5:40, 80:120])
eofs.plot.imshow(ax=ax_rgb_map)
ax_rgb_map.set_xlim(map_xlim)
ax_rgb_map.set_ylim(map_ylim)
phf.plot_juris_bounds(ax_rgb_map,
                      crs=crs,
                      strip_axes=True,
                      reset_axlims=False,
                     )
for i, row in sw_taxa.reset_index().iterrows():
    tid = row['tid']
    name = row['name']
    K = K_vals[name]
    colors = clust_colors[:K]
    ax_map = fig.add_subplot(gs[5:15, 125+(i*15):135+(i*15)])
    ax_ts = fig.add_subplot(gs[15:20, 125+(i*15):135+(i*15)])
    ax_radar = fig.add_subplot(gs[25:40, 125+(i*15):135+(i*15)],
                               projection='polar',
                              )
    # read observations
    obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
    obs = gpd.read_file(os.path.join(inat_obs_data_dir, obs_fn))
    phf.plot_flowerdate_LSP_comparison(flower_obs=obs,
                                       ax_map=ax_map,
                                       ax_ts=ax_ts,
                                       ax_radar=ax_radar,
                                       plot_crs=crs,
                                       map_xlim=map_xlim,
                                       map_ylim=map_ylim,
                                       interp_lsp_data=interp_lsp_data,
                                       neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                                       radar_alpha=0.5,
                                       radar_width_shrink_factor=0.9,
                                       colors=colors,
                                       save_scree_plot=False,
                                       name=row['name'],
                                       tid=row['tid'],
                                      )


# select and plot exemplary taxa in South African Cape
zaf_names = ['Satyrium parviflorum',
            'Pelargonium sidoides',
            'Pelargonium englerianum',
           ]
zaf_taxa = inat_mmrr_filt_adeq_n[inat_mmrr_filt_adeq_n['name'].isin(zaf_names)]
# set values of K (based on manually inspected scree plots of LSP clustering)
K_vals = {'Satyrium parviflorum': 2,
          'Pelargonium sidoides': 2,
          'Pelargonium englerianum': 2,
         }
# set map bounding box
map_xlim = [12, 36]
map_ylim = [-36, -21]
map_lim_df = pd.DataFrame({'geometry':[Point(*lims) for lims in zip(map_xlim,
                                                                    map_ylim)],
                           'idx': range(2)})
map_lim_gdf = gpd.GeoDataFrame(map_lim_df, geometry='geometry', crs=4326)
map_lim_proj_coords = map_lim_gdf.to_crs(crs).get_coordinates().values
map_xlim = map_lim_proj_coords[:, 0]
map_ylim = map_lim_proj_coords[:, 1]
ax_rgb_map = fig.add_subplot(gs[60:110, :60])
eofs.plot.imshow(ax=ax_rgb_map)
ax_rgb_map.set_xlim(map_xlim)
ax_rgb_map.set_ylim(map_ylim)
phf.plot_juris_bounds(ax_rgb_map,
                      crs=crs,
                      strip_axes=True,
                      reset_axlims=False,
                     )
for i, row in zaf_taxa.reset_index().iterrows():
    tid = row['tid']
    name = row['name']
    K = K_vals[name]
    colors = clust_colors[:K]
    ax_map = fig.add_subplot(gs[65:75, 65+(i*15):75+(i*15)])
    ax_ts = fig.add_subplot(gs[75:80, 65+(i*15):75+(i*15)])
    ax_radar = fig.add_subplot(gs[85:100, 65+(i*15):75+(i*15)],
                               projection='polar',
                              )
    # read observations
    obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
    obs = gpd.read_file(os.path.join(inat_obs_data_dir, obs_fn))
    phf.plot_flowerdate_LSP_comparison(flower_obs=obs,
                                       ax_map=ax_map,
                                       ax_ts=ax_ts,
                                       ax_radar=ax_radar,
                                       plot_crs=crs,
                                       map_xlim=map_xlim,
                                       map_ylim=map_ylim,
                                       interp_lsp_data=interp_lsp_data,
                                       neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                                       radar_alpha=0.5,
                                       radar_width_shrink_factor=0.9,
                                       colors=colors,
                                       save_scree_plot=False,
                                       name=row['name'],
                                       tid=row['tid'],
                                      )


# find and plot top candidate taxa in Eastern Brazil
# (need to be taxa with adequate and fairly balanced sampling across the
# climatic/LSP gradient described in Thomé et al. 2021 and in our paper)
countries = gpd.read_file(phf.ADM0_BOUNDS).to_crs(crs).loc[:, ['adm0_a3', 'geometry']]
subnational = gpd.read_file(phf.SELECT_ADM1_BOUNDS).to_crs(crs)
ne_brz = subnational[(subnational['NAME_1'].isin(['Maranhão',
                                                  'Tocantins',
                                                  'Piauí',
                                                  'Ceará',
                                                  'RioGrandedoNorte',
                                                  'Paraíba',
                                                  'Pernambuco',
                                                  'Alagoas',
                                                  'Sergipe',
                                                  'Bahia',
                                                 ])) &
                     (subnational['COUNTRY'] == 'Brazil')].dissolve()
se_brz = subnational[(subnational['NAME_1'].isin(['Goiás',
                                                  'DistritoFederal',
                                                  'MinasGerais',
                                                  'EspíritoSanto',
                                                  'RiodeJaneiro',
                                                  'SãoPaulo',
                                                  'Paraná',
                                                  'SantaCatarina',
                                                  'RioGrandedoSul',
                                                 ])) &
                    (subnational['COUNTRY'] == 'Brazil')].dissolve()
ne_brz_taxa = ne_brz.loc[:, ['geometry',
                             'COUNTRY']].sjoin(inat_mmrr_filt_adeq_n.to_crs(crs))
se_brz_taxa = se_brz.loc[:, ['geometry',
                             'COUNTRY']].sjoin(inat_mmrr_filt_adeq_n.to_crs(crs))
# get intersecting observation counts for each taxon with
# LSP-coefficient value <=0.05
total_obs_counts = {}
ne_brz_obs_cts = {}
se_brz_obs_cts = {}
for i, row in ne_brz_taxa.iterrows():
    if row['lsp_p']<=0.05:
        tid = row['tid']
        name = row['name']
        obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
        obs = gpd.read_file(os.path.join(inat_obs_data_dir, obs_fn))
        ct = len(ne_brz.sjoin(obs.to_crs(ne_brz.crs)))
        assert len(obs) >= ct
        print((f"\nTID {tid}: {name} has {len(obs)} total observations and "
               f"{ct} observations within the northeastern Brazil region.\n")
             )
        ne_brz_obs_cts[tid] = ct
    total_obs_counts[tid] = len(obs)
for i, row in se_brz_taxa.iterrows():
    if row['lsp_p']<=0.05:
        tid = row['tid']
        name = row['name']
        obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
        obs = gpd.read_file(os.path.join(inat_obs_data_dir, obs_fn))
        ct = len(se_brz.sjoin(obs.to_crs(se_brz.crs)))
        assert len(obs) >= ct
        print((f"\nTID {tid}: {name} has {len(obs)} total observations and "
               f"{ct} observations within the southeastern Brazil region.\n")
             )
        se_brz_obs_cts[tid] = ct
    if tid in total_obs_counts:
        assert len(obs) == total_obs_counts[tid]
    else:
        total_obs_counts[tid] = len(obs)
# get ratios or proportions of total observations
# (to find taxa with most even sampling across the region)
ne_se_brz_obs_prop_ratios = {}
for tid, total_ct in total_obs_counts.items():
    if tid in ne_brz_obs_cts and tid in se_brz_obs_cts:
        ne_prop = ne_brz_obs_cts[tid]/total_ct
        se_prop = se_brz_obs_cts[tid]/total_ct
        props = sorted([ne_prop, se_prop])
        ratio = props[0]/props[1]
        ne_se_brz_obs_prop_ratios[tid] = ratio
# get taxa that have samples across the eastern Brazil region, along with their
# sampling ratios and total sample counts in the region
e_brz_taxa = ne_brz_taxa[ne_brz_taxa['tid'].isin(ne_se_brz_obs_prop_ratios)]
e_brz_taxa.loc[:, 'reg_obs_prop_ratio'] = [ne_se_brz_obs_prop_ratios[tid] for
                                           tid in e_brz_taxa['tid'].values]
e_brz_taxa.loc[:, 'reg_obs_ct'] = [ne_brz_obs_cts[tid] + se_brz_obs_cts[tid] for
                                           tid in e_brz_taxa['tid'].values]
# subset to only those with minimum required number of samples,
# then order by sample proportion ratio and take most evenly sampled taxon
e_brz_taxa = e_brz_taxa[e_brz_taxa['reg_obs_ct'] >= min_n_inat_samps_for_demo_plots]
e_brz_taxa = e_brz_taxa.sort_values(by='reg_obs_prop_ratio', ascending=False)
# and get the SE Brazilian polygon, too
# (to clip out a handful of very distant samples)
e_brz = ne_brz.union(se_brz)

# these taxa will be plotted after landgen results, in the next section...




#------------------------------------------------------------------------------
# run landscape genetic MMRRs and visualize results
#------------------------------------------------------------------------------

# set map bounding box
map_xlim = [-70, -30]
map_ylim = [-50, 0]
map_lim_df = pd.DataFrame({'geometry':[Point(*lims) for lims in zip(map_xlim,
                                                                    map_ylim)],
                           'idx': range(2)})
map_lim_gdf = gpd.GeoDataFrame(map_lim_df, geometry='geometry', crs=4326)
map_lim_proj_coords = map_lim_gdf.to_crs(crs).get_coordinates().values
map_xlim = map_lim_proj_coords[:, 0]
map_ylim = map_lim_proj_coords[:, 1]
ax_rgb_map = fig.add_subplot(gs[125:215, :40])
eofs.plot.imshow(ax=ax_rgb_map)
ax_rgb_map.set_xlim(map_xlim)
ax_rgb_map.set_ylim(map_ylim)
phf.plot_juris_bounds(ax_rgb_map,
                      crs=crs,
                      strip_axes=True,
                      reset_axlims=False,
                     )

# analyze and plot Rhinella granulosa data
rg_l_gs_col = 45
rg_r_gs_col = 75
rg_gen_dist, rg_pts = rg.run_analysis()
ax_rg_genclust_map = fig.add_subplot(gs[130:160, rg_l_gs_col:rg_r_gs_col])
ax_rg_genclust_ts = fig.add_subplot(gs[160:170, rg_l_gs_col:rg_r_gs_col])
ax_rg_lspclust_map = fig.add_subplot(gs[175:205, rg_l_gs_col:rg_r_gs_col])
ax_rg_lspclust_ts = fig.add_subplot(gs[205:215, rg_l_gs_col:rg_r_gs_col])
phf.plot_popgen_LSP_comparison(gen_dist_mat=rg_gen_dist,
                               pts=rg_pts,
                               ax_lspclust_map=ax_rg_lspclust_map,
                               ax_lspclust_ts=ax_rg_lspclust_ts,
                               ax_genclust_map=ax_rg_genclust_map,
                               ax_genclust_ts=ax_rg_genclust_ts,
                               plot_crs=crs,
                               map_xlim=map_xlim,
                               map_ylim=map_ylim,
                               interp_lsp_data=interp_lsp_data,
                               neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                              )

# analyze and plot Xiphorhynchus fuscus data
xf_l_gs_col = 80
xf_r_gs_col = 110
xf_gen_dist, xf_pts = xf.run_analysis()
ax_xf_genclust_map = fig.add_subplot(gs[130:160, xf_l_gs_col:xf_r_gs_col])
ax_xf_genclust_ts = fig.add_subplot(gs[160:170, xf_l_gs_col:xf_r_gs_col])
ax_xf_lspclust_map = fig.add_subplot(gs[175:205, xf_l_gs_col:xf_r_gs_col])
ax_xf_lspclust_ts = fig.add_subplot(gs[205:215, xf_l_gs_col:xf_r_gs_col])
phf.plot_popgen_LSP_comparison(gen_dist_mat=xf_gen_dist,
                               pts=xf_pts,
                               ax_lspclust_map=ax_xf_lspclust_map,
                               ax_lspclust_ts=ax_xf_lspclust_ts,
                               ax_genclust_map=ax_xf_genclust_map,
                               ax_genclust_ts=ax_xf_genclust_ts,
                               plot_crs=crs,
                               map_xlim=map_xlim,
                               map_ylim=map_ylim,
                               interp_lsp_data=interp_lsp_data,
                               neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                              )

# read and combine results tables from landgen analyses, save as supp table
rg_df = pd.read_csv(os.path.join(phf.TABS_DIR,
                                 'rhinella_granulosa_MMRR_res.csv')).T
rg_df.columns = ['Rhinella granulosa']
xf_df = pd.read_csv(os.path.join(phf.TABS_DIR,
                                 'xiphorhynchus_fuscus_MMRR_res.csv')).T
xf_df.columns = ['Xiphorhynchus fuscus']
out_df = pd.concat((rg_df, xf_df), axis=1).T
out_df = out_df.loc[:, ['Intercept',
                        'Intercept(t)',
                        'Intercept(p)',
                        'LSP_dist',
                        'LSP_dist(t)',
                        'LSP_dist(p)',
                        'env_dist',
                        'env_dist(t)',
                        'env_dist(p)',
                        'geo_dist',
                        'geo_dist(t)',
                        'geo_dist(p)',
                        'R^2',
                        'F-statistic',
                        'F p-value',
                       ]
                   ]
out_df.columns = ['Intercept: coeff.',
                  'Intercept: t stat.',
                  'Intercept: P-value',
                  'LSP distance: coeff.',
                  'LSP distance: t stat.',
                  'LSP distance: P-value',
                  'Env. distance: coeff.',
                  'Env. distance: t stat.',
                  'Env. distance: P-value',
                  'Geo. distance: coeff.',
                  'Geo. distance: t stat.',
                  'Geo. distance: P-value',
                  'R^2',
                  'F stat.',
                  'F stat. P-value',
                 ]
out_df = out_df.T
out_df.to_csv(os.path.join(phf.TABS_DIR, 'TAB_SUPP_landgen_MMRR_results.csv'),
              index=True,
             )

# plot data for the two most significant iNat taxa that have substantial and well
# balanced (north:south) sampling within the same eastern Brazilian region
# (bottom right of figure, where it will appear next to the landscape genetic
# analyses for the same region)
e_brz_l_gs_cols = [110, 145]
e_brz_r_gs_cols = [125, 170]

# set values of K (based on manually inspected scree plots of LSP clustering)
K_vals = [2, # Pleroma heteromallum
          3, #Marsypianthes chamaedrys
         ]

for i, row in e_brz_taxa.reset_index().iterrows():
    K = K_vals[i]
    colors = clust_colors[:K]
    l_gs_col = e_brz_l_gs_cols[i]
    r_gs_col = e_brz_r_gs_cols[i]
    ax_map = fig.add_subplot(gs[130:160, l_gs_col:r_gs_col])
    ax_ts = fig.add_subplot(gs[160:170, l_gs_col:r_gs_col])
    ax_radar = fig.add_subplot(gs[180:210, l_gs_col:r_gs_col],
                               projection='polar',
                              )
    # read observations
    tid = row['tid']
    name = row['name']
    obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
    obs = gpd.read_file(os.path.join(inat_obs_data_dir, obs_fn))
    phf.plot_flowerdate_LSP_comparison(flower_obs=obs,
                                       ax_map=ax_map,
                                       ax_ts=ax_ts,
                                       ax_radar=ax_radar,
                                       plot_crs=crs,
                                       map_xlim=map_xlim,
                                       map_ylim=map_ylim,
                                       interp_lsp_data=interp_lsp_data,
                                       neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                                       radar_alpha=0.5,
                                       radar_width_shrink_factor=0.9,
                                       colors=colors,
                                       save_scree_plot=False,
                                       name=row['name'],
                                       tid=row['tid'],
                                      )




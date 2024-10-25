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
import matplotlib.image as mpimg
from scipy import stats
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

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/asynch/anal/phen/'))
import comparar_LSP_y_cosecha_cafetera_colombiana as cafe

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/asynch/anal/gen/rhinella/'))
import test_rhinella_granulosa as rg

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/asynch/anal/gen/xiphorhynchus/'))
import test_xiphorhynchus_fuscus as xf


#------------------------------------------------------------------------------
# set overarching params and load supporting data:
#------------------------------------------------------------------------------

# whether to interpolate missing LSP coefficients, 
# and if so, how far away can rasterio look for neighbor to interpolate from?
# coefficient raster?
interp_lsp_data = False
neigh_dist_lsp_fill_tol = 2

# plotting params for LSP time series and iNat flowering-date plots
ts_linewidth = 0.2
flower_obs_hatch_marker = '|'
flower_obs_hatch_size = 25
radar_alpha=0.5
radar_width_shrink_factor=0.9

# plotting font sizes
section_lab_fontsize=17
cbar_lab_fontsize = 16
title_fontsize = 14
taxon_fontsize = 9
axlabel_fontsize = 10

# configure figure and axes sizes
figsize = (14.5, 7)
gridspec_dims = (100, 205)
sw_cont_map_slices = (slice(2, 17),
                     slice(0, 15),
                    )
sw_photo_slices = [(slice(15, 37),
                    slice(0, 19),
                     )
                    ]
sw_scat_map_slices = [(slice(6, 36),
                       slice(20, 50),
                      ),
                     ]
sw_ts_slices = [(slice(36, 46),
                 slice(20, 50),
                ),
               ]
zaf_cont_map_slices = (slice(53, 68),
                      slice(0, 15),
                     )
zaf_photo_slices = [(slice(71, 95),
                     slice(0, 18),
                    )
                   ]
zaf_scat_map_slices = [(slice(54, 88),
                        slice(20, 50),
                       ),
                      ]
zaf_ts_slices = [(slice(88, 98),
                  slice(20, 50),
                 ),
                ]
e_brz_cont_map_slices = (slice(4, 19),
                        slice(54, 69),
                       )
e_brz_photo_slices = [(slice(19, 37),
                       slice(54, 68),
                      ),
                      (slice(58, 83),
                       slice(50, 72),
                      ),
                     ]
e_brz_genclust_map_slices = [(slice(13, 43),
                              slice(61, 91),
                             ),
                             (slice(60, 90),
                              slice(61, 91),
                             ),
                            ]
e_brz_genclust_ts_slices = [(slice(43, 53),
                             slice(61, 91),
                            ),
                            (slice(90, 100),
                             slice(61, 91),
                            ),
                            ]
e_brz_lspclust_map_slices = [(slice(13, 43),
                             slice(96, 126),
                            ),
                            (slice(60, 90),
                             slice(96, 126),
                            ),
                           ]
e_brz_lspclust_ts_slices = [(slice(43, 53),
                             slice(96, 126),
                            ),
                            (slice(90, 100),
                             slice(96, 126),
                            ),
                           ]
cafe_cont_map_slices = (slice(5, 20),
                       slice(128, 143),
                      )
cafe_photo_slices = (slice(4, 20),
                     slice(137, 167),
                    )
cafe_rgb_map_slices = (slice(20, 100),
                       slice(130, 160),
                      )
cafe_ts_slices = [(slice(13+(21*i),34+(21*i)), slice(165, 205)) for i in range(4)]

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


#------------------------------------------------------------------------------
# summarize and plot flower-phenology analysis results
#------------------------------------------------------------------------------

print(f"\n\n{'-'*80}\n\nNOW RUNNING INATURALIST ANALYSES...\n")

# create the overall figure
fig = plt.figure(figsize=figsize)
gs = fig.add_gridspec(nrows=gridspec_dims[0], ncols=gridspec_dims[1])
gs_supp = fig.add_gridspec(ncols=1, nrows=16)

# create the iNat peak-analysis supp figure
fig_supp = plt.figure(figsize=(9,15))

# plot inat flower phenology peak-analysis results
peaks_hex_filename = 'inat_hex_results.json'
inat_h3_gdf = gpd.read_file(os.path.join(inat_hex_data_dir, peaks_hex_filename))

inat_h3_gdf['prop_non1peak_signif'] = (inat_h3_gdf['prop_0peak_signif'] +
                                  inat_h3_gdf['prop_2pluspeak_signif'])
cmap='viridis'
res_cols = ['prop_0peak_signif',
            'prop_2pluspeak_signif',
            'prop_non1peak_signif',
           ]
label_dict = {'prop_non1peak_signif': 'all non-unimodal',
              'prop_0peak_signif': 'no significant peaks',
              'prop_2pluspeak_signif': '≥2 significant peaks',
             }
for i, res_col in enumerate(res_cols):
    gs_supp_imin = 5*i
    gs_supp_imax = 5*i+5+(1*(i==2))
    ax = fig_supp.add_subplot(gs_supp[gs_supp_imin:gs_supp_imax, 0])
    # add bottom axes for a colorbar
    if i == 2:
        position = 'bottom'
        pad = 0.2
        divider = make_axes_locatable(ax)
        cax = divider.append_axes(position, size='7%', pad=pad)
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
                                 vmax=1,
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
        orientation = 'horizontal'
        plt.colorbar(scalcmap, cax=cax, orientation=orientation)
        ticks = np.linspace(0, 1, 5)
        cax.set_xlabel('proportion of taxa',
                       fontdict={'fontsize': cbar_lab_fontsize},
                      )
        cax.set_xticks(ticks, ['%0.2f' % t for t in ticks], size=12)
        cax.set_ylabel('')
        cax.set_yticks(())
    # set equal aspect ratio
    ax.set_aspect('equal')

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

all_taxa_count = len(gpd.read_file('./phen/all_inat_plant_phen_taxa.csv'))
print(f"\n\n{all_taxa_count} have flowering observations in iNaturalist.")

gt50_obs_taxa_count = len(gpd.read_file('./phen/inat_flower_phen_results.json'))
print(f"\n\n{gt50_obs_taxa_count} have >= 50 flowering observations available.")

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

# correction for false discovery rate
p_vals_fdr = stats.false_discovery_control(inat_mmrr_gdf.loc[:, 'lsp_p'],
                                           method='bh')
inat_mmrr_gdf.loc[:, 'lsp_p_fdr'] = p_vals_fdr
inat_mmrr_gdf.loc[:, 'lsp_p_fdr_sig'] = p_vals_fdr <= 0.05

# sort rows by significance and then by sample size
inat_mmrr_gdf = inat_mmrr_gdf.sort_values(by=['extreme_lat_range', 'lsp_p', 'n'],
                                          ascending=[True, True, False],
                                         )

# drop extreme-latitudinal-range taxa from further analysis
inat_mmrr_filt = inat_mmrr_gdf[~inat_mmrr_gdf['extreme_lat_range']]

# save MMRR results formatted as a supplemental table
# (only taxa that are significant after FDR correction; otherwise would
# be a huge table for supps, and the full table will be in Zenodo anyhow)
supp_tab = inat_mmrr_filt[inat_mmrr_filt['lsp_p_fdr_sig']].loc[:, ['name',
                                                                   'tid',
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
                                                                   'f',
                                                                   'f_p',
                                                                   'r2',
                                                                  ]]
supp_tab.columns = ['taxon',
                    'iNat ID',
                    'n',
                    'int.',
                    'int. t',
                    'int. P',
                    'geo.',
                    'geo. t',
                    'geo. P',
                    'env.',
                    'env. t',
                    'env. P',
                    'LSP',
                    'LSP t',
                    'LSP P',
                    'F',
                    'F P',
                    'R2',
                   ]
supp_tab.to_csv(os.path.join(phf.TABS_DIR, 'TAB_SUPP_iNat_MMRR_results.csv'),
                index=False,
               )

# print & log numbers & percentages of taxa passing the different filtering & analysis stages
n_nonunimodal_taxa = np.sum(inat_res_gdf['signif_npeaks']!=1)
n_all_taxa = len(inat_res_gdf)
sig_nonunimodal_msg = (f"{np.round(100*(n_nonunimodal_taxa/n_all_taxa),1)}% "
                       "of taxa are significantly non-unimodal "
                       f"({n_nonunimodal_taxa}/{n_all_taxa}).")
npeaks_cts = inat_res_gdf.loc[:, ['signif_npeaks', 'tid']].groupby('signif_npeaks').count()
print(f"(counts of taxa by number of peaks: {npeaks_cts})")
n_taxa_w_failed_MMRR = len(inat_mmrr_df) - len(inat_mmrr_df_valid_pval)
fail_fit_MMRR_msg = (f"{n_taxa_w_failed_MMRR} non-unimodal taxa failed "
                      "to fit a valid MMRR model.")
n_taxa_w_broad_distrs = len(inat_mmrr_gdf) - len(inat_mmrr_filt)
too_broad_distr_msg = (f"{n_taxa_w_broad_distrs} non-unimodal taxa dropped "
                        "because of broad, equator-crossing distributions.")
n_MMRR_taxa = len(inat_mmrr_filt)
n_MMRR_taxa_LSP_signif_0p05 = np.sum(inat_mmrr_filt['lsp_p']<=0.05)
n_MMRR_taxa_LSP_signif_bonf_corr = np.sum(inat_mmrr_filt['lsp_p_sig_bonf_corr'])
n_MMRR_taxa_LSP_signif_fdr_corr = np.sum(inat_mmrr_filt['lsp_p_fdr_sig'])
sig_MMRR_LSP_msg = (f"{np.round(100*(n_MMRR_taxa_LSP_signif_0p05/n_MMRR_taxa),1)}% "
                     "of non-unimodal taxa for which MMRRs were fitted "
                     "had significant LSP-distance coefficients (P<=0.05) "
                    f"({n_MMRR_taxa_LSP_signif_0p05}/{n_MMRR_taxa}).")
sig_MMRR_LSP_msg_bonf_corr = (
                f"{np.round(100*(n_MMRR_taxa_LSP_signif_bonf_corr/n_MMRR_taxa),1)}% "
                 "of non-unimodal taxa for which MMRRs were fitted "
                 "had significant LSP-distance coefficients after "
                f"Bonferroni correction (P-value <={np.round(p_bonf_corr, 5)}) "
                f"({n_MMRR_taxa_LSP_signif_bonf_corr}/{n_MMRR_taxa}).")
sig_MMRR_LSP_msg_fdr_corr = (
                f"{np.round(100*(n_MMRR_taxa_LSP_signif_fdr_corr/n_MMRR_taxa),1)}% "
                 "of non-unimodal taxa for which MMRRs were fitted "
                 "had significant LSP-distance coefficients after "
                f"FDR correction (P-value <=0.05) "
                f"({n_MMRR_taxa_LSP_signif_fdr_corr}/{n_MMRR_taxa}).")

print(f"{fail_fit_hull_msg}\n\n")
print(f"{sig_nonunimodal_msg}\n\n")
print(f"{fail_fit_MMRR_msg}\n\n")
print(f"{too_broad_distr_msg}\n\n")
print(f"{sig_MMRR_LSP_msg}\n\n")
print(f"{sig_MMRR_LSP_msg_bonf_corr}\n\n")
print(f"{sig_MMRR_LSP_msg_fdr_corr}\n\n")

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
        f.write(f"{sig_MMRR_LSP_msg_fdr_corr}\n\n")

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
        new_cols['prop_lsp_signif'].append(
            np.mean(stats.false_discovery_control(gdf_intsxn['lsp_p'],
                                                  method='bh')<=0.05))
        new_cols['prop_lsp_signif_0p05'].append(np.mean(gdf_intsxn['lsp_p']<=0.05))
        new_cols['mean_geo_p'].append(np.mean(gdf_intsxn['geo_p']))
        new_cols['mean_geo'].append(np.mean(gdf_intsxn['geo']))
        new_cols['prop_geo_signif'].append(
            np.mean(stats.false_discovery_control(gdf_intsxn['geo_p'],
                                                  method='bh')<=0.05))
        new_cols['mean_int_p'].append(np.mean(gdf_intsxn['int_p']))
        new_cols['mean_int'].append(np.mean(gdf_intsxn['int']))
        new_cols['prop_int_signif'].append(
            np.mean(stats.false_discovery_control(gdf_intsxn['int_p'],
                                                  method='bh')<=0.05))
        new_cols['mean_f_p'].append(np.mean(gdf_intsxn['f_p']))
        new_cols['prop_mod_signif'].append(
            np.mean(stats.false_discovery_control(gdf_intsxn['f_p'],
                                                  method='bh')<=0.05))
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
# plot demonstrative iNaturalist taxa
#------------------------------------------------------------------------------

# drop small sample sizes (i.e., where too many LSP pixels were masked out)
inat_mmrr_filt_adeq_n = inat_mmrr_filt[inat_mmrr_filt['n'] >=min_n_inat_samps_for_demo_plots]

# set cluster colors
flower_asynch_colors = np.array(['#c4ba78', # tan
                                 '#2cab87', # turquoise
                                ])

def plot_taxon_photo(ax, name):
    """
    plot a photo of the taxon, taken from
    the available CC BY and CC BY-NC photos on iNaturalist
    """
    filename = name.replace(' ', '_') + '.png'
    if name == 'Xiphorhynchus fuscus':
        filename = os.path.splitext(filename)[0] + '_horizontal.png'
    filepath = os.path.join(phf.EXTERNAL_INAT_DATA_DIR,
                            'photos',
                            filename,
                           )
    img = mpimg.imread(filepath)
    # make white pixels (i.e. pixels with 1s in the first 3 layers) transparent
    img[:,:,3] *= (img[:,:,:3].sum(axis=2)!=3)
    # swap x and y axes, if necessary
    if 'horizontal' in filename:
        img = img.swapaxes(0, 1)
    ax.imshow(img)
    ax.set_aspect('equal')
    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_axis_off()
    ax.patch.set_alpha(0)


def plot_map_bbox(map_xlims,
                  map_ylims,
                  ax,
                  box_buff_pct=0.05,
                  alpha=1.0,
                  linewidth=0.5,
                  label=None,
                  labelxpad=0.1,
                  labelypad=0.1,
                  size=12,
                 ):
    """
    adds the bounding box indicated by the given x- and y-lims to the given
    axes
    """
    if box_buff_pct != 0:
        x_buff = box_buff_pct * (np.abs(np.diff(map_xlims)))
        y_buff = box_buff_pct * (np.abs(np.diff(map_ylims)))
        buff = np.mean((x_buff, y_buff))
        map_xlims = [map_xlims[0] - buff, map_xlims[1] + buff]
        map_ylims = [map_ylims[0] - buff, map_ylims[1] + buff]
    xs = np.array(map_xlims)[[0, 0, 1, 1, 0]]
    ys = np.array(map_ylims)[[0, 1, 1, 0, 0]]
    ax.plot(xs, ys, '-k', alpha=alpha, linewidth=linewidth, clip_on=False)
    if label is not None:
        assert labelxpad is not None
        assert labelypad is not None
        ax.text(map_xlims[0] + labelxpad * (np.abs(np.diff(map_xlims))),
                map_ylims[0] + labelypad * (np.abs(np.diff(map_ylims))),
                label,
                weight='bold',
                size=size,
               )


def get_projected_maplims(map_xlims, map_ylims, crs):
    '''
    return the map x- and y-limits, provided in unproject lat-lon
    as equivalent coordinates in the indicated crs
    '''
    map_lims_df = pd.DataFrame({'geometry':[Point(*lims) for lims in zip(map_xlims,
                                                                        map_ylims)],
                               'idx': range(2)})
    map_lims_gdf = gpd.GeoDataFrame(map_lims_df, geometry='geometry', crs=4326)
    map_lims_proj_coords = map_lims_gdf.to_crs(crs).get_coordinates().values
    map_xlims_proj = map_lims_proj_coords[:, 0]
    map_ylims_proj = map_lims_proj_coords[:, 1]
    return map_xlims_proj, map_ylims_proj


def plot_continental_reference_map(ax,
                                   map_xlims,
                                   map_ylims,
                                   box_xlims,
                                   box_ylims,
                                   crs=crs,
                                  ):
    '''
    plot a continental-scale map of black jurisdictional bounds,
    as a reference map for a regional analysis
    '''
    ax.set_xlim(map_xlims)
    ax.set_ylim(map_ylims)
    ax.set_aspect('equal')
    phf.plot_juris_bounds(ax,
                          crs=crs,
                          strip_axes=True,
                          reset_axlims=False,
                          lev0_linewidth=0.2,
                          lev1_linewidth=0.05,
                         )
    ax.set_axis_off()
    plot_map_bbox(box_xlims,
                  box_ylims,
                  ax=ax,
                  label=None,
                  linewidth=1,
                 )


def plot_focal_inat_taxa(mmrr_res_gdf,
                         taxa,
                         ax_cont_map,
                         scatter_map_axs,
                         ts_axs,
                         photo_axs,
                         flow_obs_axs,
                         map_xlims,
                         map_ylims,
                         cont_map_xlims,
                         cont_map_ylims,
                         flow_obs_plot_type='stack',
                         radar_alpha=0.5,
                         radar_width_shrink_factor=0.9,
                         save_scree_plot=False,
                         set_title=False,
                        ):
    """
    function for identically visualizing a series of example taxa
    """
    assert len(taxa) == len(scatter_map_axs) == len(ts_axs)
    assert flow_obs_axs is None or len(taxa) == len(flow_obs_axs)
    # reproject the given map lims to the required CRS
    (map_xlims_proj,
     map_ylims_proj) = get_projected_maplims(map_xlims, map_ylims, crs)
    (cont_map_xlims_proj,
     cont_map_ylims_proj) = get_projected_maplims(cont_map_xlims,
                                                  cont_map_ylims,
                                                  crs,
                                                 )
    if ax_cont_map is not None:
        plot_continental_reference_map(ax_cont_map,
                                       cont_map_xlims_proj,
                                       cont_map_ylims_proj,
                                       map_xlims_proj,
                                       map_ylims_proj,
                                       crs=crs,
                                      )
    ct = 0
    for taxon, K in taxa.items():
        # plot the photo first (so its axes don't overlap anything else)
        plot_taxon_photo(photo_axs[ct], taxon)

        tax_dict = mmrr_res_gdf[mmrr_res_gdf['name'] == taxon].iloc[0,:]
        tid = tax_dict['tid']
        name = tax_dict['name']
        print(f"\n\tplotting {name}...")
        ax_map = scatter_map_axs[ct]
        ax_ts = ts_axs[ct]
        if flow_obs_axs is not None:
            ax_flow_obs = flow_obs_axs[ct]
        else:
            ax_flow_obs = None
        # read observations
        obs_fn = f"TID_{tid}_{name.replace(' ', '_')}.json"
        obs = gpd.read_file(os.path.join(inat_obs_data_dir, obs_fn))
        # drop points outside max boundaries
        len_b4 = len(obs)
        in_bounds = []
        for pt in obs.to_crs(4326).get_coordinates().values:
            in_bounds.append((map_xlims[0]<=pt[0]<=map_xlims[1] and
                              map_ylims[0]<=pt[1]<=map_ylims[1]))
        obs = obs[in_bounds]
        len_af = len(obs)
        print(f"{len_b4-len_af} points dropped outside map bounds")
        phf.plot_flowerdate_LSP_comparison(flower_obs=obs,
                                           ax_map=ax_map,
                                           ax_ts=ax_ts,
                                           ax_flow_obs=ax_flow_obs,
                                           K=K,
                                           plot_crs=crs,
                                           map_xlim=map_xlims_proj,
                                           map_ylim=map_ylims_proj,
                                           interp_lsp_data=interp_lsp_data,
                                           neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                                           ts_linewidth=ts_linewidth,
                                           flower_obs_plot_type=None,
                                           flower_obs_hatch_marker=flower_obs_hatch_marker,
                                           flower_obs_hatch_size=flower_obs_hatch_size,
                                           radar_alpha=radar_alpha,
                                           radar_width_shrink_factor=radar_width_shrink_factor,
                                           colors=flower_asynch_colors,
                                           save_scree_plot=save_scree_plot,
                                           name=name,
                                           tid=tid,
                                          )
        # set axis titles after all other components have been plotted
        if set_title:
            ax_map.set_title("$" + name.replace(' ', '\ ') + "$",
                             fontdict={'fontsize': title_fontsize},
                            )
        else:
            ax_map.set_title('')
        ct+=1


# . . . . . . . . . . . . . . . . . . . .
# plot example taxon in SW USA and Mexico
# . . . . . . . . . . . . . . . . . . . .
# set taxa and their K values (based on manual inspection of scree plots)
sw_taxa_clust_Ks = {'Xanthisma spinulosum': 2}
# set up axes objects and params for plotting
map_xlims = [-115, -105]
map_ylims = [24.66, 38]
cont_map_xlims = [-135, -60]
cont_map_ylims = [17, 52]
ax_cont_map = fig.add_subplot(gs[sw_cont_map_slices[0],
                                sw_cont_map_slices[1]])
scatter_map_axs = [fig.add_subplot(gs[sw_scat_map_slices[i][0],
                                      sw_scat_map_slices[i][1]]) for i in range(len(sw_taxa_clust_Ks))]
ts_axs = [fig.add_subplot(gs[sw_ts_slices[i][0],
                             sw_ts_slices[i][1]]) for i in range(len(sw_taxa_clust_Ks))]
photo_axs = [fig.add_subplot(gs[sw_photo_slices[i][0],
                                sw_photo_slices[i][1]]) for i in range(len(sw_taxa_clust_Ks))]
plot_focal_inat_taxa(mmrr_res_gdf=inat_mmrr_filt_adeq_n,
                     taxa=sw_taxa_clust_Ks,
                     ax_cont_map=ax_cont_map,
                     scatter_map_axs=scatter_map_axs,
                     ts_axs=ts_axs,
                     photo_axs=photo_axs,
                     flow_obs_axs=None,
                     map_xlims=map_xlims,
                     map_ylims=map_ylims,
                     cont_map_xlims=cont_map_xlims,
                     cont_map_ylims=cont_map_ylims,
                     radar_alpha=radar_alpha,
                     radar_width_shrink_factor=radar_width_shrink_factor,
                     save_scree_plot=False,
                    )
[ax.set_ylabel('scaled LSP',
               labelpad=7,
               fontdict={'fontsize': axlabel_fontsize},
             ) for ax in ts_axs]
[ax.set_aspect('equal') for ax in scatter_map_axs]


# . . . . . . . . . . . . . . . .
# plot example taxon in S. Africa
# . . . . . . . . . . . . . . . .
# set taxa and their K values (based on manual inspection of scree plots)
zaf_taxa_clust_Ks = {'Satyrium parviflorum': 2}
# set map bounding box
map_xlims = [17, 31]
map_ylims = [-35.5, -22]
cont_map_xlims = [-20, 59]
cont_map_ylims = [-39, 39]
ax_cont_map = fig.add_subplot(gs[zaf_cont_map_slices[0],
                                zaf_cont_map_slices[1]])
scatter_map_axs = [fig.add_subplot(gs[zaf_scat_map_slices[i][0],
                                      zaf_scat_map_slices[i][1]]) for i in range(len(zaf_taxa_clust_Ks))]
ts_axs = [fig.add_subplot(gs[zaf_ts_slices[i][0],
                             zaf_ts_slices[i][1]]) for i in range(len(zaf_taxa_clust_Ks))]
photo_axs = [fig.add_subplot(gs[zaf_photo_slices[i][0],
                                zaf_photo_slices[i][1]]) for i in range(len(sw_taxa_clust_Ks))]
plot_focal_inat_taxa(mmrr_res_gdf=inat_mmrr_filt_adeq_n,
                     taxa=zaf_taxa_clust_Ks,
                     ax_cont_map=ax_cont_map,
                     scatter_map_axs=scatter_map_axs,
                     ts_axs=ts_axs,
                     photo_axs=photo_axs,
                     flow_obs_axs=None,
                     map_xlims=map_xlims,
                     map_ylims=map_ylims,
                     cont_map_xlims=cont_map_xlims,
                     cont_map_ylims=cont_map_ylims,
                     radar_alpha=radar_alpha,
                     radar_width_shrink_factor=radar_width_shrink_factor,
                     save_scree_plot=False,
                    )
[ax.set_ylabel('scaled LSP',
               labelpad=7,
               fontdict={'fontsize': axlabel_fontsize},
             ) for ax in ts_axs]
[ax.set_aspect('equal') for ax in scatter_map_axs]


#---------------------------------------------
# run coffee-harvest analysis and plot results
#---------------------------------------------

print(f"\n\n{'-'*80}\n\nNOW RUNNING COFFEA ANALYSIS...\n")

# set up axgs
# add black box to bigger regional map to indicate focal region
map_xlims = [-7.45e6, -6.863e6]
map_ylims = [1.38e5, 1.46e6]
SAm_cont_map_xlims = [-7.816e6, -3.180e6]
SAm_cont_map_ylims = [-6.815e6, 1.828e6]
ts_axs = [fig.add_subplot(gs[cafe_ts_slices[i][0],
                             cafe_ts_slices[i][1]]) for i in range(len(cafe_ts_slices))]
rgb_map_ax = fig.add_subplot(gs[cafe_rgb_map_slices[0],
                                cafe_rgb_map_slices[1]])
# run and plot the analysis
n_perms = 1000
cafe.run_analysis(rgb_map_ax,
                  ts_axs,
                  map_xlims,
                  map_ylims,
                  n_perms=n_perms,
                  eofs_alpha=0.75,
                  region_sample_marker_size=40,
                 )
# add map of bigger region
ax_cont_map = fig.add_subplot(gs[cafe_cont_map_slices[0],
                                cafe_cont_map_slices[1]])
plot_continental_reference_map(ax_cont_map,
                               SAm_cont_map_xlims,
                               SAm_cont_map_ylims,
                               map_xlims,
                               map_ylims,
                               crs=crs,
                              )
# add photo of species
ax_photo = fig.add_subplot(gs[cafe_photo_slices[0],
                              cafe_photo_slices[1],
                             ])
plot_taxon_photo(ax_photo, 'Coffea arabica')


#------------------------------------------------------------------------------
# run landscape genetic MMRRs and visualize results
#------------------------------------------------------------------------------

landgen_colors=np.array(['#648fff', # blue
                         '#e84138', # red
                        ])

# set map bounding box
map_xlims = [-48, -34]
map_ylims = [-24.3, -2.8]
map_lims_df = pd.DataFrame({'geometry':[Point(*lims) for
                                lims in zip(map_xlims, map_ylims)],
                            'idx': range(2)})
map_lims_gdf = gpd.GeoDataFrame(map_lims_df, geometry='geometry', crs=4326)
map_lims_proj_coords = map_lims_gdf.to_crs(crs).get_coordinates().values
map_xlims = map_lims_proj_coords[:, 0]
map_ylims = map_lims_proj_coords[:, 1]
ax_cont_map = fig.add_subplot(gs[e_brz_cont_map_slices[0],
                                e_brz_cont_map_slices[1]])
plot_continental_reference_map(ax_cont_map,
                               SAm_cont_map_xlims,
                               SAm_cont_map_ylims,
                               map_xlims,
                               map_ylims,
                               crs=crs,
                              )


print(f"\n\n{'-'*80}\n\nNOW RUNNING RHINELLA LANDGEN ANALYSIS...\n")

# analyze and plot Rhinella granulosa data
rg_gen_dist, rg_pts = rg.run_analysis()
ax_rg_genclust_map = fig.add_subplot(gs[e_brz_genclust_map_slices[0][0],
                                        e_brz_genclust_map_slices[0][1]])
ax_rg_genclust_ts = fig.add_subplot(gs[e_brz_genclust_ts_slices[0][0],
                                       e_brz_genclust_ts_slices[0][1]])
ax_rg_lspclust_map = fig.add_subplot(gs[e_brz_lspclust_map_slices[0][0],
                                        e_brz_lspclust_map_slices[0][1]])
ax_rg_lspclust_ts = fig.add_subplot(gs[e_brz_lspclust_ts_slices[0][0],
                                       e_brz_lspclust_ts_slices[0][1]])
phf.plot_popgen_LSP_comparison(gen_dist_mat=rg_gen_dist,
                               pts=rg_pts,
                               ax_lspclust_map=ax_rg_lspclust_map,
                               ax_lspclust_ts=ax_rg_lspclust_ts,
                               ax_genclust_map=ax_rg_genclust_map,
                               ax_genclust_ts=ax_rg_genclust_ts,
                               K=2,
                               colors=landgen_colors,
                               plot_crs=crs,
                               map_xlim=map_xlims,
                               map_ylim=map_ylims,
                               interp_lsp_data=interp_lsp_data,
                               neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                              )
# add photo
ax_rg_photo = fig.add_subplot(gs[e_brz_photo_slices[0][0],
                                 e_brz_photo_slices[0][1]])
plot_taxon_photo(ax_rg_photo, 'Rhinella granulosa')
# set title and labels last 
ax_rg_genclust_map.set_title('')
ax_rg_genclust_ts.set_ylabel('scaled LSP',
                             labelpad=7,
                              fontdict={'fontsize': axlabel_fontsize},
                             )


print(f"\n\n{'-'*80}\n\nNOW RUNNING XIPHORHYNCHUS LANDGEN ANALYSIS...\n")

# analyze and plot Xiphorhynchus fuscus data
xf_gen_dist, xf_pts = xf.run_analysis()
ax_xf_genclust_map = fig.add_subplot(gs[e_brz_genclust_map_slices[1][0],
                                        e_brz_genclust_map_slices[1][1]])
ax_xf_genclust_ts = fig.add_subplot(gs[e_brz_genclust_ts_slices[1][0],
                                       e_brz_genclust_ts_slices[1][1]])
ax_xf_lspclust_map = fig.add_subplot(gs[e_brz_lspclust_map_slices[1][0],
                                        e_brz_lspclust_map_slices[1][1]])
ax_xf_lspclust_ts = fig.add_subplot(gs[e_brz_lspclust_ts_slices[1][0],
                                       e_brz_lspclust_ts_slices[1][1]])
phf.plot_popgen_LSP_comparison(gen_dist_mat=xf_gen_dist,
                               pts=xf_pts,
                               ax_lspclust_map=ax_xf_lspclust_map,
                               ax_lspclust_ts=ax_xf_lspclust_ts,
                               ax_genclust_map=ax_xf_genclust_map,
                               ax_genclust_ts=ax_xf_genclust_ts,
                               K=2,
                               colors=landgen_colors,
                               plot_crs=crs,
                               map_xlim=map_xlims,
                               map_ylim=map_ylims,
                               interp_lsp_data=interp_lsp_data,
                               neigh_dist_lsp_fill_tol=neigh_dist_lsp_fill_tol,
                              )
# add photo
ax_xf_photo = fig.add_subplot(gs[e_brz_photo_slices[1][0],
                                 e_brz_photo_slices[1][1]])
plot_taxon_photo(ax_xf_photo, 'Xiphorhynchus fuscus')
# set title and labels last
ax_xf_genclust_map.set_title('')
ax_xf_genclust_ts.set_ylabel('scaled LSP',
                             labelpad=7,
                              fontdict={'fontsize': axlabel_fontsize},
                             )
# force all map axes to equal aspect ratios
for ax in [ax_rg_genclust_map, ax_rg_lspclust_map,
           ax_xf_genclust_map, ax_xf_lspclust_map]:
    ax.set_aspect('equal')


# add a single, giant transparent axes over top of everything,
# then use that to add section labels and dividers, etc
ax_meta = fig.add_subplot(gs[:, :])
ax_meta.patch.set_alpha(0)
ax_meta.axis('off')
ax_meta.set_xticks(())
ax_meta.set_yticks(())
ax_meta.set_xlabel('')
ax_meta.set_ylabel('')
ax_meta.set_title('')
ax_meta.text(0.0075,
             0.98,
             'A. Flowering asynchrony',
             weight='bold',
             size=section_lab_fontsize,
             clip_on=False,
            )
ax_meta.text(0.265,
             0.98,
             'B. Genetic isolation by asynchrony',
             weight='bold',
             size=section_lab_fontsize,
             clip_on=False,
            )
ax_meta.text(0.635,
             0.98,
             'C. Harvest asynchrony',
             weight='bold',
             size=section_lab_fontsize,
             clip_on=False,
            )
ax_meta.text(0.37,
             0.89,
             'genetic\nclusters',
             ha='center',
             rotation=0,
             size=title_fontsize,
            )
ax_meta.text(0.541,
             0.89,
             'LSP\nclusters',
             ha='center',
             rotation=0,
             size=title_fontsize,
            )
ax_meta.text(0.04,
             0.595,
             [*sw_taxa_clust_Ks][0].replace(' ', '\n'),
             ha='center',
             fontdict={'fontsize': taxon_fontsize,
                       'style': 'italic',
                      },
            )
ax_meta.text(0.04,
             0.005,
             [*zaf_taxa_clust_Ks][0].replace(' ', '\n'),
             ha='center',
             fontdict={'fontsize': taxon_fontsize,
                       'style': 'italic',
                      },
             clip_on=False,
            )
ax_meta.text(0.30,
             0.63,
             'Rhinella\ngranulosa',
             ha='center',
             fontdict={'fontsize': taxon_fontsize,
                       'style': 'italic',
                      },
            )
ax_meta.text(0.30,
             0.13,
             'Xiphorhynchus\nfuscus',
             ha='center',
             fontdict={'fontsize': taxon_fontsize,
                       'style': 'italic',
                      },
            )
ax_meta.text(0.74,
             0.79,
             'Coffea\narabica',
             ha='center',
             fontdict={'fontsize': taxon_fontsize,
                       'style': 'italic',
                      },
            )
ax_meta.plot([0.255, 0.255],
             [-0.2, 1.2],
             linewidth=0.7,
             color='black',
             alpha=0.7,
             clip_on=False,
             zorder=0,
            )
ax_meta.plot([0.6275, 0.6275],
             [-0.2, 1.2],
             linewidth=0.7,
             color='black',
             alpha=0.7,
             clip_on=False,
             zorder=0,
            )
ax_meta.set_xlim(0, 1)
ax_meta.set_ylim(0, 1)

# adjust subplots and save
fig.subplots_adjust(hspace=0,
                    wspace=0,
                    left=0,
                    right=0.99,
                    bottom=0.04,
                    top=0.98,
                   )
fig.savefig(os.path.join(phf.FIGS_DIR, 'FIG_inat_landgen_coffea_results.png'),
            dpi=600,
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


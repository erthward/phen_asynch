#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
from sklearn.cluster import KMeans
from shapely.geometry import Point
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf


def run_analysis(rgb_map_ax,
                 ts_axs,
                 map_xlims,
                 map_ylims,
                 eofs_alpha=0.75,
                 region_sample_marker_size=40,
                ):
    # set seed, for reproducibility of permutation testing
    np.random.seed(0)

    #--------------------------------------------------------------------------
    # data prep

    # load digitized cosecha points
    data_dir = os.path.join(phf.REPO_DIR, 'src/asynch/anal/phen/')
    patt = '(VERDE)|(AMARILLO)|(ANARANJADO)|(MORADO)'
    files = [f for f in os.listdir(data_dir) if re.search(patt, f)]
    # reorder files so that they plot in the overal N-S order of the clusters
    NS_color_order = ['VERDE', 'ANARANJADO', 'MORADO', 'AMARILLO']
    order = [np.argwhere([re.search(c, f) for f in files]) for c in NS_color_order]
    order = np.stack(order).ravel()
    files = [*np.array(files)[order]]
    # I've named files based on the colors used in the figure I digitized, but we
    # should instead choose colors that better highlight the overarching N-S
    # shift in LSP shape and timing (but that are still distinct enough to clearly
    # highlight the very short-distance regional discontinuities)
    cmap = mpl.cm.viridis
    rgbas = cmap(np.linspace(0, 1, 4))
    hexes = [mpl.colors.rgb2hex(v, keep_alpha=True) for v in rgbas]
    colors = {'VERDE': hexes[3],
              'ANARANJADO': hexes[2],
              'MORADO': hexes[1],
              'AMARILLO': hexes[0],
             }
    gdfs = []
    for f in files:
        df = pd.read_csv(os.path.join(data_dir, f), header=None)
        df.loc[:, 'geometry'] = [Point(*row.values) for i, row in df.iterrows()]
        principal = re.search('(?<=principal)[A-Z]{4}', f).group()
        mitaca = re.search('(?<=mitaca)[A-Z]{2}', f)
        df.loc[:, 'principal'] = principal
        if mitaca is not None:
            df.loc[:, 'mitaca'] = mitaca.group()
        else:
            df.loc[:, 'mitaca'] = np.nan
        color = re.search(patt, f).group()
        df.loc[:, 'color'] = colors[color]
        gdf = gpd.GeoDataFrame(df, geometry='geometry', crs=4326)
        gdfs.append(gdf)
    gdf = pd.concat(gdfs)
    # reset index, to facilitate dropping points with missing LSP data
    gdf = gdf.reset_index()

        # NOTE: numerical doy range associated with the major and minor seasons,
    #       adjusted down by 1 to account for January 1st being doy 0
    season_xlims = {'MAMJ': [59, 180],
                    'SOND': [243, 364],
                    'AM': [90, 150],
                    'ON': [273, 333],
                   }


    #--------------------------------------------------------------------------
    # time series extraction and plotting

    all_keep = []
    all_lsp = []
    for i, color in enumerate(colors.values()):
        ax = ts_axs[i]
        # extract LSP ts at each point
        subgdf = gdf.loc[gdf['color'] == color]
        pts = subgdf.get_coordinates().values
        lsp_ts = phf.get_raster_info_points(phf.COEFFS_FILE,
                                            pts,
                                            'ts',
                                            standardize=True,
                                            fill_nans=False,
                                            fill_tol=None,
                                           )
        assert np.all(lsp_ts.shape == (len(pts), 365))
        # drop NA points
        keep = np.sum(pd.isnull(lsp_ts), axis=1) == 0
        all_keep.extend([*subgdf.index.values[keep]])
        print((f"\n\t{np.round(100*np.mean(np.invert(keep)), 1)}% of points "
               f"in cluster {i} are missing LSP data...\n"))
        lsp_ts = lsp_ts[keep, :]
        all_lsp.append(lsp_ts)
        mid = np.median(lsp_ts, axis=0)
        # plot regions' low and high percentiles
        pctile_lo = np.percentile(lsp_ts, 10, axis=0)
        pctile_hi = np.percentile(lsp_ts, 90, axis=0)
        pctile_poly_coords = np.vstack((np.array([*zip(range(365), pctile_lo)]),
                            np.array([*zip([*range(365)][::-1], pctile_hi[::-1])])))
        pctile_poly = Polygon(pctile_poly_coords)
        patchcoll = PatchCollection([pctile_poly],
                                    alpha=0.3,
                                    color=color,
                                    edgecolor='k',
                                   )
        ax.add_collection(patchcoll)
        # plot the median as a darker line on top
        ax.plot(range(len(mid)),
                mid,
                linestyle='-',
                color=color,
                alpha=1,
                zorder=1,
                linewidth=3,
               )
        principal = np.unique(subgdf['principal'])
        assert len(principal) == 1
        principal = principal[0]
        ax.plot(season_xlims[principal],
                [-2.3]*2,
                color=color,
                linewidth=5,
                zorder=0,
               )
        if np.sum(pd.isnull(subgdf['mitaca'])) == 0:
            mitaca = np.unique(subgdf['mitaca'])
            assert len(mitaca) == 1
            mitaca = mitaca[0]
            ax.plot(season_xlims[mitaca],
                    [-2.3]*2,
                    color=color,
                    linewidth=5,
                    zorder=0,
                   )

        ax.set_xlim(0, 364)
        ax.set_ylim(-2.5, 2)
        ax.set_yticks(())
        if i == 2:
            ax.set_ylabel(f"{' '*27}scaled LSP", fontdict={'fontsize': 10})
        else:
            ax.set_ylabel('')
        if i == 3:
            ax.set_xticks([0, 90, 181, 273, 364],
                          ['Jan', 'Apr', 'Jul', 'Oct', 'Jan'],
                         )
            ax.tick_params(labelsize=8)
        else:
            ax.set_xlabel('')
            ax.set_xticks(())


    #--------------------------------------------------------------------------
    # time series clustering and statistical anaylsis of cluster agreement

    # keep only points with LSP data
    gdf = gdf.loc[all_keep, :]

        # concatenate all LSP time series into a single array
    all_lsp_arr = np.vstack(all_lsp)
    # cluster the LSP time series at all the sampling points
    clusts = KMeans(n_clusters=4).fit(all_lsp_arr)
    # calculate binary matrices indicating whether each pair of points has the same
    # or different labels, both for the officially mapped coffee-harvest regions
    # and for our clustering labels
    print("\n\tnow calculating Jaccard index...\n")
    reg_same = (np.ones([len(gdf)]*2) * np.nan).astype(np.int16)
    clust_same = (np.ones([all_lsp_arr.shape[0]]*2) * np.nan).astype(np.int16)
    assert np.all(reg_same.shape == clust_same.shape)
    for i in range(len(gdf)):
        for j in range(len(gdf)):
            if i == j:
                reg_same[i,j] = 1
                clust_same[i, j] = 1
            else:
                reg_same[i,j] = reg_same[j, i] = (gdf['color'].iloc[i] ==
                                                  gdf['color'].iloc[j])
                clust_same[i,j] = clust_same[j, i] = (clusts.labels_[i] ==
                                                      clusts.labels_[j])
    # extract the upper triangular values indicating pairwise label agreements
    reg_same_vals = reg_same[np.triu_indices(reg_same.shape[0], k=1)]
    clust_same_vals = clust_same[np.triu_indices(clust_same.shape[0], k=1)]
    # check length is correct
    expected_len = ((reg_same.shape[0] **2) - reg_same.shape[0])/2
    assert len(reg_same_vals) == len(clust_same_vals) == expected_len
    # get vectors indicating whether point-pairs have the same label in both
    # datasets, the same label only in the coffee harvest regions, or the same
    # label only in our clustering results
    same_both = reg_same_vals & clust_same_vals
    same_reg_only = reg_same_vals & ~clust_same_vals
    same_clust_only = ~reg_same_vals & clust_same_vals
    sum_same_either = np.sum((same_both, same_reg_only, same_clust_only))
    # check that the sum of all three of those vectors equals the sum of all
    # point-pairs that are the same in at least one of the two datasets
    # (i.e, either the same in the coffee harvest regions or the same in our
    # clustering results)
    assert sum_same_either == np.sum(reg_same_vals | clust_same_vals)
    # use those vectors to calculate the Jaccard index for our clustering vis-a-vis
    # the coffee harvest regions' sample points
    jaccard_ind = np.sum(same_both)/sum_same_either
    assert 0 <= jaccard_ind <= 1
    # also use the same approach to simulate 999 Jaccard index values from datasets
    # in which the clustering is run after the assignment of LSP time series to
    # points has been permuted
    perm_jaccard_inds = []
    n_perms = 1000
    print("\n\tnow calculating permutated Jaccard indices...\n")
    for i in range(n_perms):
        if (i+1) % 100 == 0:
            print(f"\n\t{np.round(100*((i+1)/n_perms),1)}% complete")
        perm_inds = np.random.choice(a=[*range(all_lsp_arr.shape[0])],
                                     size=all_lsp_arr.shape[0],
                                     replace=False,
                                    )
        all_lsp_arr_perm = all_lsp_arr[perm_inds, :]
        assert np.all(all_lsp_arr_perm.shape == all_lsp_arr.shape)
        perm_clusts = KMeans(n_clusters=4).fit(all_lsp_arr_perm)
        perm_clust_same = (np.ones([all_lsp_arr_perm.shape[0]]*2) * np.nan).astype(np.int16)
        assert np.all(reg_same.shape == perm_clust_same.shape)
        for i in range(len(gdf)):
            for j in range(len(gdf)):
                if i == j:
                    perm_clust_same[i, j] = 1
                else:
                    perm_clust_same[i,j] = perm_clust_same[j, i] = (perm_clusts.labels_[i] ==
                                                                    perm_clusts.labels_[j])
        perm_clust_same_vals = perm_clust_same[np.triu_indices(perm_clust_same.shape[0], k=1)]
        expected_len = ((reg_same.shape[0] **2) - reg_same.shape[0])/2
        assert len(perm_clust_same_vals) == expected_len
        perm_same_both = reg_same_vals & perm_clust_same_vals
        perm_same_reg_only = reg_same_vals & ~perm_clust_same_vals
        perm_same_clust_only = ~reg_same_vals & perm_clust_same_vals
        perm_sum_same_either = np.sum((perm_same_both, perm_same_reg_only, perm_same_clust_only))
        assert perm_sum_same_either == np.sum(reg_same_vals | perm_clust_same_vals)
        perm_jaccard_ind = np.sum(perm_same_both)/perm_sum_same_either
        perm_jaccard_inds.append(perm_jaccard_ind)
    # use results to calculate the empirical P-value of our observed Jaccard index
    assert len(perm_jaccard_inds) == n_perms
    assert np.all([0 <= pji <= 1 for pji in perm_jaccard_inds])
    Pval = np.mean(jaccard_ind <= np.array(perm_jaccard_inds))

    print(("\n\nJaccard index between coffee-harvest region assignment and "
           f"K-means clustering (k=4): {np.round(jaccard_ind, 3)} "
           f"(permutation-based P-value={str(np.round(Pval, 3))};"
           f"i.e., << {1/n_perms}) "
           "(range of permuted Jaccard values: "
           f"[{np.round(np.min(perm_jaccard_inds), 3)}, "
           f"{np.round(np.max(perm_jaccard_inds), 3)}])\n\n"))


    #--------------------------------------------------------------------------
    # mapping of harvest region sampling points and LSP diversity

    # read in the unfolded EOFs file
    # (NOTE: avoids the 'color smudging' artefact caused by folding the EOFS across
    # the ITCZ, which isn't an issue in other focal/regional maps but is confounding here)
    eofs = rxr.open_rasterio(phf.EOFS_FILE)[:3].rio.set_crs(4326)
    for i in range(3):
        eofs[i] = (eofs[i] - np.nanmin(eofs[i]))/(np.nanmax(eofs[i])-np.nanmin(eofs[i]))
    eofs = eofs.rio.reproject(8857)
    gdf_plot = gdf.to_crs(8857)


    # map regions, as concave hulls, on top of EOFs map
    rgb_map_ax.set_aspect('equal')
    eofs.plot.imshow(ax=rgb_map_ax,
                     alpha=eofs_alpha,
                    )
    for color in np.unique(gdf['color']):
        subgdf = gdf_plot[gdf_plot['color'] == color]
        subgdf.plot(ax=rgb_map_ax,
                    color=color,
                    edgecolor='k',
                    marker='.',
                    markersize=region_sample_marker_size,
                    linewidth=0.4,
                    zorder=2,
                   )

    phf.plot_juris_bounds(rgb_map_ax,
                          lev0_linewidth=0.3,
                          lev0_zorder=1,
                          lev1=False,
                          crs=8857,
                          strip_axes=True,
                          reset_axlims=False,
                         )
    rgb_map_ax.set_xlim(*map_xlims)
    rgb_map_ax.set_ylim(*map_ylims)


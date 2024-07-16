"""
test our fitted LSP patterns against all
GenBank Xiphorhynchus fuscus data used
in the only eastern Brazilian analysis in
Qunitero et al. 2014 (DOI 10.1086/677261),
originally from Cabanne et al. 2008
(DOI 10.1016/j.ympev.2008.09.013)
by running an MMRR and using our LSP
patterns to test the Asynchrony of Seasons Hypothesis
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import pyproj
import rasterio as rio
import rioxarray as rxr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import re
import pprint

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf
from MMRR import MMRR


# set the numpy.random seed
np.random.seed(1)

# whether to interpolate missing LSP coefficients, 
# and if so, how far away can rasterio look for neighbor to interpolate from?
# coefficient raster?
interp_lsp_data = False
neigh_dist_lsp_fill_tol = 5

# how many MMRR permutations to use
MMRR_nperm = 999

def run_analysis():
    # read genetic distance matrix
    print('\nreading genetic data...')
    gen_dist_filepath = os.path.join(phf.REPO_DIR,
                                     'src/asynch/anal/gen/xiphorhynchus/',
                                     'all_xiphorhynchus_fuscus_dist_IAN.csv',
                                    )
    gen_dist = pd.read_csv(gen_dist_filepath)
    # set first column as index
    gen_dist.set_index(gen_dist.columns[0], inplace=True)
    # reformat row and col names to reflect just the voucher code
    gen_dist.index = [i.split('.')[0] for i in gen_dist.index]
    gen_dist.columns = [c.split('.')[0] for c in gen_dist.columns]
    # replace nullified upper half of matrix with symmetric values
    gen_dist.values[np.triu_indices(len(gen_dist))] = gen_dist.T.values[np.triu_indices(len(gen_dist))]
    # and replace nullified self-distances with zeros
    for i in range(len(gen_dist)):
        gen_dist.iloc[i, i] = 0
    # check voucher IDs on rows and columns are unique and identical
    assert len(np.unique(gen_dist.index)) == len(gen_dist)
    assert len(np.unique(gen_dist.columns)) == len(gen_dist.columns)
    assert np.all(np.array([*gen_dist.index]) == np.array([*gen_dist.columns]))
    # check diagonal is all zeros
    assert np.all([gen_dist.iloc[i, i] == 0 for i in range(gen_dist.shape[0])])
    # check matrix is symmetric
    assert np.all(gen_dist.T == gen_dist)
    # read CSV with sample voucher codes and locations
    print('\nreading geographic data...')
    geo_filepath = os.path.join(phf.REPO_DIR,
                                'src/asynch/anal/gen/xiphorhynchus/',
                                'all_xiphorhynchus_fuscus_locs.csv',
                               )
    geo = pd.read_csv(geo_filepath)
    # calculate decimal degree columns
    # NOTE: only missing vals are geo coord secs in most rows, so replace w/ 0s
    geo = geo.replace(np.nan, 0)
    # NOTE: multiply all by -1 because all in S and W hemispheres
    geo['lat'] = -1*(geo['lat_d'] + ((geo['lat_m'] + (geo['lat_s']/60))/60))
    geo['lon'] = -1*(geo['lon_d'] + ((geo['lon_m'] + (geo['lon_s']/60))/60))
    # check we have matching voucher IDs in both datasets
    assert set(geo['id'].values) == set(gen_dist.index)

    # set up the geographic, environmental, and LSP distance matrices
    geo_dist = np.ones([geo.shape[0]]*2) * np.nan

    # loop over voucher codes in gen_dist, to ensure that rows & cols
    # match the rows & cols in gen_dist
    print('\ncalculating geographic distances...')
    g = pyproj.Geod(ellps='WGS84')
    pts = []
    for i, vouch_i in enumerate(gen_dist.index):
        lon_i, lat_i = geo.loc[geo['id'] == vouch_i, ['lon', 'lat']].values[0]
        pts.append([lon_i, lat_i])
        for j, vouch_j in enumerate(gen_dist.columns):
            # calculate geographic distance
            if i == j:
                dist_ij = 0
            else:
                lon_j, lat_j = geo.loc[geo['id'] == vouch_j, ['lon', 'lat']].values[0]
                (az_ij, az_ji, dist_ij) = g.inv(lon_i, lat_i, lon_j, lat_j)
            geo_dist[i, j] = dist_ij
            geo_dist[j, i] = dist_ij
    pts = np.array(pts)

    # calculate environmental distance across all Worldclim vals
    print('\ncalculating environmental distances...')
    nodata_val = rio.open(phf.BIOCLIM_INFILEPATHS[0]).nodata
    env_dist = phf.calc_pw_clim_dist_mat(pts,
                                         nodata_val=nodata_val,
                                        )

    # calculate LSP distance
    # NOTE: needs to be based on standardized time series,
    #       to focus on timing rather than fitted NIRV values
    print('\ncalculating LSP distances...')
    lsp_dist = phf.get_raster_info_points(phf.COEFFS_STRICT_FILE,
                                          pts,
                                          'ts_pdm',
                                          standardize=True,
                                          fill_nans=interp_lsp_data,
                                          fill_tol=neigh_dist_lsp_fill_tol,
                                         )
    print((f"\n\n{np.round(100*np.mean(pd.isnull(lsp_dist)), 2)} % "
            "of lsp_dist consists of missing values"))

    # check all diagonals are zeros
    assert np.all([geo_dist[i, i] == 0 for i in range(geo_dist.shape[0])])
    assert np.all([env_dist[i, i] == 0 for i in range(env_dist.shape[0])])
    assert np.all([lsp_dist[i, i] == 0 for i in range(lsp_dist.shape[0])])

    # drop rows and cols in gen_dist, geo_dist, and env_dist
    # if they feature NaNs in lsp_dist
    # (happens because our raster has NaNs where pixels were filtered out)
    print('\ndropping rows with missing LSP data...')
    missing_lsp = (np.sum(pd.isnull(lsp_dist), axis=0) == lsp_dist.shape[0] - 1)
    missing_lsp_vouch = gen_dist.index[np.where(missing_lsp)].values
    print(("\n\nThe following voucher IDs occur at locations where LSP "
           "data is missing and was not interpolated:\t"
           f"{'   '.join(missing_lsp_vouch)}\n\n"))

    fig, ax = plt.subplots(1,1)
    phf.plot_juris_bounds(ax,
                          lev1=False,
                          lev0_linewidth=1,
                          lev0_alpha=1,
                          crs=4326,
                          strip_axes=False,
                         )
    missing_lsp_color = ['blue' if v not in missing_lsp_vouch else 'red' for
                                                        v in geo['id'].values]
    geo['missing_lsp'] = missing_lsp_color
    ax.scatter(x=geo['lon'], y=geo['lat'], c=geo['missing_lsp'], s=25, alpha=0.5)
    ax.set_title(('red points missing data in our LSP dataset'
                  f'\n({np.round(100*np.mean(missing_lsp), 2)}% of all sample sites)')
                )
    ax.set_xlim(-57.5, -34)
    ax.set_ylim(-30, 0)
    fig.savefig(os.path.join(phf.FIGS_DIR,
                             'xiphorhynchus_fuscus_sites_missing_LSP_data.png'), dpi=300)

    gen_dist = gen_dist.values
    gen_dist = gen_dist[~missing_lsp].T[~missing_lsp]
    geo_dist = geo_dist[~missing_lsp].T[~missing_lsp]
    env_dist = env_dist[~missing_lsp].T[~missing_lsp]
    lsp_dist = lsp_dist[~missing_lsp].T[~missing_lsp]
    # check symmetry and shape of all matrices
    for i, mat in enumerate([gen_dist, geo_dist, env_dist, lsp_dist]):
        assert np.all(mat == mat.T), f"matrix {i} failed!"
        assert np.all(mat.shape == gen_dist.shape), f"matrix {i} failed!"


    # run the MMRR model and print results
    print('\nrunning MMRR model...')
    res = MMRR(Y=gen_dist,
               X=[geo_dist, env_dist, lsp_dist],
               Xnames=['geo_dist', 'env_dist', 'LSP_dist'],
               # NOTE: MMRR will standardize lower-triangular distance values, and thus
               #       returns coefficient values as beta-coefficients
               standardize=True,
               intercept=True,
               nperm=MMRR_nperm,
              )
    print('\nMMRR results:\n')
    pprint.pprint(res)

    # save MMRR results to a table
    res_df = pd.DataFrame({k:[v] for k, v in res.items()})
    res_df.to_csv(os.path.join(phf.TABS_DIR,
                               'xiphorhynchus_fuscus_MMRR_res.csv'),
                  index=False,
                 )

    # return genetic distance matrix, and points that aren't missing LSP data,
    # for plotting of results
    pts_for_plot = pts[~missing_lsp, :]
    return gen_dist, pts_for_plot

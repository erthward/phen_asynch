"""
test our fitted LSP patterns against the data from
Thomé et al. 2021 (DOI 10.1038/s41437-021-00460-7)
by redoing their MMRR test but using our LSP
patterns to test the Asynchrony of Seasons hypothesis
in lieu of their monthly climatic LDA variables
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
import re
import sys
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
    # NOTE: simple Euclidean genetic distance matrix calculated using adegenet in R
    gen_dist_filepath = os.path.join(phf.REPO_DIR,
                                     'src/asynch/anal/gen/rhinella/dryad_archive/',
                                     '2-Mantel_test/Rgranulosa_Mantel.csv',
                                    )
    gen_dist = pd.read_csv(gen_dist_filepath)
    # reformat row and col names to reflect just the voucher code
    # (without the lingering haplotype integer)
    gen_dist.index = [i.split('_')[0] for i in gen_dist.index]
    gen_dist.columns = [c.split('_')[0] for c in gen_dist.columns]
    # check voucher IDs on rows and columns are unique and identical
    assert len(np.unique(gen_dist.index)) == len(gen_dist)
    assert len(np.unique(gen_dist.columns)) == len(gen_dist.columns)
    assert np.all(np.array([*gen_dist.index]) == np.array([*gen_dist.columns]))
    # check diagonal is all zeros
    assert np.all([gen_dist.iloc[i, i] == 0 for i in range(gen_dist.shape[0])])

    # read CSV with sample voucher codes and locations
    print('\nreading geographic data...')
    geo_filepath = os.path.join(phf.REPO_DIR,
                                'src/asynch/anal/gen/rhinella/',
                                '41437_2021_460_MOESM1_ESM_TABLES1.csv',
                               )
    geo = pd.read_csv(geo_filepath)
    # get rid of white space in column names and values
    geo.columns = [c.strip() for c in geo.columns]
    # rename lon and lat
    geo.rename(mapper={'Longitude': 'lon', 'Latitude': 'lat'}, axis=1, inplace=True)
    geo['voucher'] = [v.strip() for v in geo['voucher'].values]
    # NOTE: MTCT1468 in gen_dist rows/columns must match
    #       MTCT1368 in geo['voucher'], as they're the only
    #       ones that reciprocally have no match in each other.
    #       so fix that!
    geo.loc[geo['voucher'] == 'MTCT1368', 'voucher'] = 'MTCT1468'
    # check we have matching voucher IDs in both datasets
    assert set(geo['voucher'].values) == set(gen_dist.index)


    # set up the geographic, environmental, and seasonal distance matrices
    geo_dist = np.ones([geo.shape[0]]*2) * np.nan

    # loop over voucher codes, ordered by gen_dist rows,
    # to ensure that rows & cols match the rows & cols in gen_dist
    vouchers = np.array([*gen_dist.index])
    print('\ncalculating geographic distances...')
    g = pyproj.Geod(ellps='WGS84')
    pts = []
    for i, vouch_i in enumerate(vouchers):
        lon_i, lat_i = geo.loc[geo['voucher'] == vouch_i, ['lon', 'lat']].values[0]
        pts.append([lon_i, lat_i])
        for j, vouch_j in enumerate(gen_dist.columns):
            # calculate geographic distance
            if i == j:
                dist_ij = 0
            else:
                lon_j, lat_j = geo.loc[geo['voucher'] == vouch_j, ['lon', 'lat']].values[0]
                (az_ij, az_ji, dist_ij) = g.inv(lon_i, lat_i, lon_j, lat_j)
            geo_dist[i, j] = dist_ij
            geo_dist[j, i] = dist_ij


    pts = np.array(pts)
    # calculate environmental distance for just the 4 variables they list
    # (BIO1, BIO4, BIO12, BIO15)
    print('\ncalculating environmental distances...')
    nodata_val = rio.open(phf.BIOCLIM_INFILEPATHS[0]).nodata
    env_dist = phf.calc_pw_clim_dist_mat(pts,
                                         nodata_val=nodata_val,
                                         vars_list=['bio_1',
                                                    'bio_4',
                                                    'bio_12',
                                                    'bio_15',
                                                   ],
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
    missing_lsp_vouch = vouchers[np.where(missing_lsp)]
    print(("\n\nThe following voucher IDs occur at locations where LSP "
           "data is missing and was not interpolated:\t"
           f"{'   '.join(missing_lsp_vouch)}\n\n"))
    fig, ax = plt.subplots(1,1)
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    world.plot(color='none', ax=ax)
    missing_lsp_color = ['blue' if v not in missing_lsp_vouch else 'red' for
                                                        v in geo['voucher'].values]
    geo['missing_lsp'] = missing_lsp_color
    ax.scatter(x=geo['lon'], y=geo['lat'], c=geo['missing_lsp'], s=25, alpha=0.5)
    ax.set_title(('red points missing data in our LSP dataset'
                  f'\n({np.round(100*np.mean(missing_lsp), 2)}% of all sample sites)')
                )
    ax.set_xlim(-57.5, -34)
    ax.set_ylim(-30, 0)
    fig.savefig(os.path.join(phf.FIGS_DIR,
                             'rhinella_granulosa_sites_missing_LSP_data.png'), dpi=300)

    gen_dist = gen_dist.values
    gen_dist = gen_dist[~missing_lsp].T[~missing_lsp]
    geo_dist = geo_dist[~missing_lsp].T[~missing_lsp]
    env_dist = env_dist[~missing_lsp].T[~missing_lsp]
    lsp_dist = lsp_dist[~missing_lsp].T[~missing_lsp]
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
                               'rhinella_granulosa_MMRR_res.csv'),
                  index=False,
                 )

    # return genetic distance matrix, and points that aren't missing LSP data,
    # for plotting of results
    pts_for_plot = pts[~missing_lsp, :]
    return gen_dist, pts_for_plot


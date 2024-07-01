#!/usr/bin/env python
# coding: utf-8
# calc_LSP_EOFs.py

import numpy as np
import pandas as pd
import geopandas as gpd
from rasterio.plot import reshape_as_raster
import xarray as xr
import rioxarray as rxr
import os
import sys
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from eofs.xarray import Eof

sys.path.insert(1, '/global/home/users/drewhart/seasonality/seasonal_asynchrony/etc/')
import phen_helper_fns as phf


###################
# BEHAVIORAL PARAMS
###################
# save results?
save_res = True

# dataset params
dataset = 'NIRv'
masking_mode = 'default'
mask_filename_ext_dict = {'strict': '_STRICTMASK',
                          'default': ''}
mask_filename_ext = mask_filename_ext_dict[masking_mode]

# plotting params
subplots_adj_left=0.05
subplots_adj_bottom=0.1
subplots_adj_right=0.95
subplots_adj_top=0.9
subplots_adj_wspace=0.2
subplots_adj_hspace=0.4

# min and max x and y values, to optionally subset analysis to a region
min_x = None
max_x = None
min_y = None
max_y = None


# data dir on laptop
if os.getcwd().split('/')[1] == 'home':
    data_dir = os.path.join('/media/deth/SLAB/diss/3-phn/final_maps_and_results/',
                            mask_filename_ext)
# data dir on savio
else:
    data_dir = os.path.join('/global/scratch/users/drewhart/seasonality/GEE_outputs',
                            dataset + mask_filename_ext)

# set the seed
set_seed = True
if set_seed:
    seed = 1
    np.random.seed(seed)

# standardize each ts to itself?
# NOTE: it makes sense to do this, since I'm only interested in timing,
#       and otherwise (especially based on Alex Turner's results in CA) I
#       expect the first EOF will largely reflect global (i.e.,
#       cross-study-area) variation in magnitude of fitted values;
#       nonetheless, setting a flag for this so that I can check that
#       expectation and check sensitivity to this decision
standardize_ts = True

# define the standardization function
def standardize(arr):
    return (arr-np.mean(arr))/np.std(arr)


# latitude weights to use?
#lat_weights = None
#lat_weights = 'cos'
lat_weights = 'sqrt_cos'

# number of top EOFs to use?
neofs = 4

# pctiles and counts of example time series to plot?
ts_top_pctile = 95
ts_bot_pctile = 95
n_ts_to_plot = 1000


#################################
# GET ARRAY OF FITTED TIME SERIES
#################################
# read global NIRv coeffs
coeffs = rxr.open_rasterio(os.path.join(data_dir, '%s%s_coeffs.tif' % (dataset, mask_filename_ext)))

# subset global raster to study area (if all Nones then not subsetted!)
# NOTE: max_y and min_y flipped because y res negative in CRS transform
coeffs = coeffs.sel(x=slice(min_x, max_x), y=slice(max_y, min_y))

# get coords arrays
X, Y = np.meshgrid(coeffs.x, coeffs.y)

# filepath for the 365 x X x Y array of fitted time series
ts_arr_filepath = os.path.join(data_dir, f'fitted_ts_array_{dataset}{mask_filename_ext}.npy')

# if that file doesn't exist then create it
if not os.path.isfile(ts_arr_filepath):
    # create empty time-series array for EOF analysis
    ts_arr = np.ones((365, coeffs.shape[1], coeffs.shape[2]), dtype=np.float32) * np.nan

    # make the harmonic regression's design matrix
    dm = phf.make_design_matrix()

    # get i and j values for non-null pixels that need time series calculated
    I, J = np.where(pd.notnull(coeffs[0,:,:]))

    # get the time series for non-null pixel
    # NOTE: coeffs.shape[0] == 5, one band for each regression coeff
    for n in range(len(I)):
        i = I[n]
        j = J[n]
        coeffs_vec = coeffs[:, i, j].values
        ts = np.sum(coeffs_vec * dm, axis=1)
        # standardize time series [0,1], if desired
        # NOTE: if not, pretty certain that first EOF will largely capture
        #       global (i.e., across full subsetted extent) variation
        #       in fitted magnitude
        if standardize_ts and not np.any(np.isnan(ts)):
            ts = standardize(ts).flatten()
        assert ts.shape == (365,)
        ts_arr[:, i, j] = ts
        if n%10000 == 0:
            print(f"\n\n\t{np.round(100*((n+1)/len(I)), 1)}% complete")
    # once complete, save this to a simple numpy array file
    # (for now, anyhow; might be worth saving as a big geospatial file eventually?)
    if save_res:
        print('\n\n\tNOW SAVING TIME SERIES CUBE...')
        np.save(ts_arr_filepath, ts_arr)
# otherwise read it in
else:
    # load saved data and reshape it from 2d to 3d
    ts_arr = np.load(ts_arr_filepath)
    assert np.all(ts_arr.shape == (365, coeffs.shape[1], coeffs.shape[2]))


#########
# RUN EOF
#########
# calculate weights array requested
if lat_weights == 'cos':
    weights = np.cos(np.deg2rad(Y))
    weights /= weights.mean()
elif lat_weights == 'sqrt_cos':
    weights = np.sqrt(np.cos(np.deg2rad(Y)))
    weights /= weights.mean()
else:
    weights = None

# coerce ts array to rio xarray obj
ts_da = xr.DataArray([coeffs[0,:,:]*np.nan]*365)
ts_da.attrs = coeffs.attrs
ts_da.attrs['long_name'] = ['d%i' % i for i in range(1, 366)]
ts_da = ts_da.rename({'dim_0': 'time',
                      'dim_1': 'y',
                      'dim_2': 'x',
                      })
ts_da = ts_da.assign_coords({'time': range(1, 366),
                             'y': coeffs.y.values,
                             'x': coeffs.x.values,
                             })
ts_da = ts_da.rio.write_crs(4326)
ts_da.rio.set_crs(4326)
ts_da.loc[:,:,:] = ts_arr

# use empirical orthogonal functions to collapse global ts into
# main modes of variation
solver = Eof(ts_da, weights=weights)

# grab the first n EOFs
eofs = solver.eofsAsCorrelation(neofs=neofs)

# grab the PCs
pcs = solver.pcs(npcs=neofs, pcscaling=1)

# grab percent variances of EOFs
var_pcts = solver.varianceFraction(neofs)

# reconstruct the field using just the selected top EOFs
ts_recon = solver.reconstructedField(neofs)

# write eofs to file, if requested
if save_res:
    tif_filename = '%s_%i_EOFs_%s%s%s.tif' % (dataset,
                                              neofs,
                                              lat_weights + 'wts',
                                              '_standts' * standardize_ts,
                                              mask_filename_ext,
                                             )
    eof_res_for_file = eofs.to_dataset('mode')
    eof_res_for_file = eof_res_for_file.rename_vars(
                                    {i: 'eof%i' % i for i in range(neofs)})
    eof_res_for_file.rio.to_raster(os.path.join(data_dir, tif_filename),
                       dtype=np.float32,
                       tags={'eof%i_pctvar' % i:
                             str(var_pcts.values[i]) for i in range(neofs)},
                      )

# save the PC values and percents variance explained to a CSV, for use in supplemental figure
pc_dict = {f"EOF{neof}_{float(np.round(var_pcts[neof]*100, 2))}pct": pcs.sel(mode=neof).values for neof in range(4)}
pc_df = pd.DataFrame.from_dict(pc_dict)
pc_df_filepath = os.path.join(data_dir, f'{dataset}{mask_filename_ext}_EOF_PC_table.csv')
pc_df.to_csv(pc_df_filepath, index=False)

# save a simple scree plot, for later reference 
scree_fig = plt.figure(figsize=(4,4))
ax = scree_fig.add_subplot(111)
ax.plot(solver.varianceFraction(10)*100)
ax.set_title('pct variance explained by first 10 EOFs')
ax.set_xlabel('EOF number')
ax.set_ylabel('pct variance explained')
scree_fig_filepath = os.path.join(data_dir, f'{dataset}{mask_filename_ext}_scree_fig.png')
scree_fig.savefig(scree_fig_filepath, dpi=500)

#!/usr/bin/env python
# coding: utf-8
# calc_global_LSP_EOFs_and_clusts.py

import numpy as np
import pandas as pd
import geopandas as gpd
from rasterio.plot import reshape_as_raster
import xarray as xr
import rioxarray as rxr
import os
import sys
import pickle
from collections import deque
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from eofs.xarray import Eof
from sklearn.cluster import MiniBatchKMeans

sys.path.insert(1,
                '/global/home/users/drewhart/seasonality/seasonal_asynchrony/src/etc/')
import phen_helper_fns as phf


###########################################
# BEHAVIORAL PARAMS AND SUPPORTING DATASETS
###########################################
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

# data dir on laptop
if os.getcwd().split('/')[1] == 'home':
    data_dir = os.path.join('/media/deth/SLAB/diss/3-phn/final_maps_and_results/',
                            mask_filename_ext)
# data dir on savio
else:
    data_dir = os.path.join('/global/scratch/users/drewhart/seasonality/GEE_outputs_and_derivatives',
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

# number of clusters to produce globally?
# NOTE: chosen by scree plot inspection
k_global = 9

# load ITCZ shapefile
# NOTE: digitized from Li and Zeng 2005, as reproduced in Zhisheng et al. 2015
itcz = gpd.read_file('/global/scratch/users/drewhart/seasonality/itcz/ITCZ_li_zeng_2005_digitized.shp')


#################################
# GET ARRAY OF FITTED TIME SERIES
#################################
# read global NIRv coeffs
coeffs = rxr.open_rasterio(os.path.join(data_dir, '%s%s_coeffs.tif' % (dataset, mask_filename_ext)))

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
scree_fig_filepath = os.path.join(data_dir, f'{dataset}{mask_filename_ext}_EOF_scree_plot.png')
scree_fig.savefig(scree_fig_filepath, dpi=500)


########################
# RUN K-MEANS CLUSTERING
########################

def run_kmeans_clust(data,
                     K,
                     batch_size=40,
                     n_init=10,
                     max_no_improvement=10,
                    ):
    '''
    run mini-batch K-means clustering and return key results
    '''
    clust = MiniBatchKMeans(init="k-means++",
                            n_clusters=K,
                            batch_size=batch_size,
                            n_init=n_init,
                            max_no_improvement=max_no_improvement,
                            verbose=0,
                           )
    clust.fit(data)
    inertia = clust.inertia_
    centers = clust.cluster_centers_
    labels = clust.labels_
    return labels, centers, inertia


k_max=13
batch_size=40
n_init=10
max_no_improvement=10
seed=1

# rotate all time series by 182 days if they're below the mean ITCZ
mean_itcz = itcz.loc[itcz['time'] == 'mean']
linestring = mean_itcz['geometry'].iloc[0]
itcz_lons = [coord[0] for coord in linestring.coords]
itcz_lats = [coord[1] for coord in linestring.coords]

# get raster dimensions' coordinates
rast_lons = coeffs.x.values
rast_lats = coeffs.y.values

# for each raster column (i.e., longitude)
for j, rast_lon in enumerate(rast_lons):
    # find the ITCZ point closest to it in longitude
    itcz_pt_idx = np.argmin(np.abs(itcz_lons-rast_lon))
    # get that ITCZ point's latitude
    itcz_lat = itcz_lats[itcz_pt_idx]
    # find the raster row closest to that in latitude
    rast_row_idx = np.argmin(np.abs(rast_lats-itcz_lat))
    # get all row indices below that
    inds_to_rotate = [*range(rast_row_idx+1, coeffs.shape[1])]
    # do rotation
    for i in inds_to_rotate:
        ts_raw = ts_arr[:, i, j]
        ts_deque = deque(ts_raw)
        ts_deque.rotate(182)
        ts_arr[:, i, j] = np.array(ts_deque)

# get fitted annual time series for all non-null pixels 
ts_arr_melt = ts_arr.reshape(365,np.prod(ts_arr.shape[1:])).T
ts_arr_melt_nonnan_inds = np.where(pd.notnull(ts_arr_melt[:, 0]))[0]
ts_arr_melt_nonnan = ts_arr_melt[ts_arr_melt_nonnan_inds, :]
assert ts_arr_melt_nonnan.shape == (len(ts_arr_melt_nonnan_inds), 365)
# run clustering, creating scree plot and saving results for the chosen
# value of k
all_labels = {}
all_centers = {}
inertias = []
for k in range(1, k_max):
    print(f"\n\tclustering global LSP time series into {k} clusters...\n")
    np.random.seed(1)
    # NOTE: using mini-batch K-means because it's faster
    labs, cents, inertia = run_kmeans_clust(data=ts_arr_melt_nonnan,
                                            K=k,
                                            batch_size=batch_size,
                                            n_init=n_init,
                                            max_no_improvement=max_no_improvement,
                                           )
    inertias.append(inertia)
    all_labels[k] = labs
    all_centers[k] = cents

labels = all_labels[k_global]
centers = all_centers[k_global]

fig_scree, ax = plt.subplots(1)
ax.plot(range(1, k_max), inertias)
ax.set_xlabel('n clusters')
ax.set_ylabel('inertia')
ax.set_title('global', fontdict={'fontsize': 24})
fig_scree.savefig(os.path.join(data_dir,
        '%s%s_LSP_clust_scree_GLOBAL.png' % (dataset,
                                             mask_filename_ext)), dpi=300)
# compose and save array of cluster assignments
cluster_arr_melt = ts_arr_melt[:, 0]*np.nan
cluster_arr_melt[ts_arr_melt_nonnan_inds] = labels
assert len([v for v in np.unique(cluster_arr_melt) if pd.notnull(v)]) == k_global
cluster_arr = cluster_arr_melt.T.reshape(ts_arr.shape[1:])
cluster_arr_path = os.path.join(data_dir,
    '%s%s_LSP_%i_clusts_arr_GLOBAL.txt' % (dataset, mask_filename_ext, k_global))
np.savetxt(cluster_arr_path, cluster_arr)
# pickle the cluster centers
centers_filepath = (os.path.join(data_dir,
    '%s%s_LSP_%i_clusts_cents_GLOBAL.pkl' % (dataset, mask_filename_ext, k_global)))
with open(centers_filepath, 'wb') as f:
    pickle.dump(centers, f)


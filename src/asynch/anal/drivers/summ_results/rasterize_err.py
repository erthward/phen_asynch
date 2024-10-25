import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
from affine import Affine
from collections import Counter as C
import os

# set working directory

# get path for savio...
if os.getcwd().split('/')[1] == 'global':
    data_dir = "/global/scratch/users/drewhart/seasonality/asynch_drivers_analysis/"
# ... or on laptop
else:
    data_dir = "/media/deth/SLAB/diss/3-phn/final_maps_and_results/rf/"

# loop over vars and neigh_rads (in km)
for var in ['NIRv', 'SIF']:
    for neigh_rad in [50, 100, 150]:
        print('\n\nPROCESSING VAR %s, NEIGH RAD %i...\n\n' % (var, neigh_rad))
        for coords_as_covars in ['y', 'n']:

            # get CSV of preds and errors
            csv = pd.read_csv(os.path.join(data_dir,
                'rf_full_preds_%sCOORDS_%s_%ikm.csv' % (coords_as_covars,
                                                        var, neigh_rad)))

            # get CSV with unprojected coordinates, to match them up to
            csv_coords = pd.read_csv(os.path.join(data_dir,
                'asynch_model_all_vars_prepped_%s_%ikm.csv' % (var, neigh_rad)))

            # convert lon, lat to floats
            csv.x = csv_coords.x
            csv.y = csv_coords.y

            # clear coords CSV from memory
            del csv_coords

            # get sorted unique lon and lat values
            xs = np.sort(np.unique(csv.x))
            #xs = np.sort([*set(csv.x)])
            ys = np.sort(np.unique(csv.y))
            #ys = np.sort([*set(csv.y)])

            # get stepwise diffs between sorted lon and lat values
            x_diffs = C(np.round(xs[1:] - xs[:-1], 3))
            y_diffs = C(np.round(ys[1:] - ys[:-1], 3))

            # get minimum diffs
            x_res = np.min([diff for diff in [*x_diffs.keys()] if diff > 0])
            y_res = np.min([diff for diff in [*y_diffs.keys()] if diff > 0])
            # NOTE: MAKE YRES NEGATIVE!
            y_res *= -1

            # get min and max lon and lat values
            x_min = np.min(np.round(xs, 3))
            x_max = np.max(np.round(xs, 3))
            # NOTE: REVERSE Y_MIN AND Y_MAX BC OF NEGATIVE RES
            y_min = np.max(np.round(ys, 3))
            y_max = np.min(np.round(ys, 3))

            # create affine transform array from those values
            transform = Affine(x_res, 0.0, x_min, 0.0, y_res, y_min)

            # make meshgrid for output raster
            X, Y = np.meshgrid(np.arange(x_min, x_max+x_res, x_res),
                               np.arange(y_min, y_max+y_res, y_res))

            # make error raster
            col = 'err'
            print('\n' + '='*80 + '\nNOW RASTERIZING %s:\n\n' % col)
            # make placeholder array for raster values
            vals = np.ones(X.shape) * np.nan
            for n, row in csv.iterrows():
                if n % 100000 == 0:
                    print('\nrow number %i of %i...\n\n' % (n, len(csv)))
                i = np.where(np.round(Y[:, 0], 3) == np.round(row['y'], 3))[0][0]
                j = np.where(np.round(X[0, :], 3) == np.round(row['x'], 3))[0][0]
                vals[i,j] = row[col]
            # cast as DataArray object
            da = xr.DataArray(vals,
                              coords = {'y': Y[:,0], 'x':X[0,:]},
                              dims=['y', 'x'],
                              attrs={'var': col})
            # convert to raster object
            da = da.rio.write_crs(4326, inplace=True)
            da = da.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
            da = da.rio.write_coordinate_system(inplace=True)
            da = da.rio.write_transform(transform, inplace=True)

            # write to file
            da.rio.to_raster(os.path.join(data_dir,
                'err_map_%sCOORDS_%s_%ikm.tif' % (coords_as_covars,
                                                  var, neigh_rad)))

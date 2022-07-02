import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
from affine import Affine
from collections import Counter as C
import os

# set working directory
data_dir = '/media/deth/SLAB/seasonality/results/'

# get CSV Shapley results
csv = pd.read_csv(os.path.join(data_dir, 'rf_SHAP_vals_w_coords.csv'))

# get sorted unique lon and lat values
xs = np.sort(np.unique(csv.x))
ys = np.sort(np.unique(csv.y))

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

# make all SHAP-value rasters
for col in [col for col in csv.columns if col not in ['x', 'y']]:
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
    da.rio.to_raster(os.path.join(data_dir, 'SHAP_map_' + col + '.tif'))

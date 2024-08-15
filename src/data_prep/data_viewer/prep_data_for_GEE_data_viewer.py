#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf

# input SRS (for datasets other than the prepped EOFs in 8857)
# (unprojected lat,lon)
in_crs = 4326

# output SRS
# (Web Mercator, since it plays well with GEE)
out_crs = 3857

# load coefficients
coeffs = rxr.open_rasterio(phf.COEFFS_FILE, masked=True)
assert coeffs.rio.crs == in_crs
coeffs = coeffs.rio.reproject(out_crs)
assert coeffs.rio.crs == out_crs

# load the raw EOFS
eofs = rxr.open_rasterio(phf.EOFS_FILE, masked=True)
eofs = eofs.rio.set_crs(in_crs)
assert eofs.rio.crs == in_crs
eofs = eofs.rio.reproject_match(coeffs)
assert eofs.rio.crs == out_crs
# replace the nodata values with np.nan
eofs = eofs.where(eofs < 1e36, np.nan).rio.set_nodata(np.nan)

# load the scaled, ITCZ-folded EOFs, for data viz of phen significant MMRR results
eofs_prepped = rxr.open_rasterio(phf.EOFS_PREPPED_FILE)
eofs_prepped = eofs_prepped.rio.write_crs(8857)
eofs_prepped = eofs_prepped.rio.reproject_match(eofs)
assert eofs_prepped.rio.crs == out_crs
# replace the nodata values with np.nan
eofs_prepped = eofs_prepped.where(eofs_prepped < 1e36, np.nan).rio.set_nodata(np.nan)

# load asynchrony
neigh = 100
asynch = rxr.open_rasterio(phf.ASYNCH_FILES[neigh], masked=True)
assert asynch.rio.crs == in_crs
asynch = asynch.rio.reproject_match(eofs_prepped)
assert asynch.rio.crs == out_crs

# stack them all and write to file
# NOTE: EOFs, asynch, then coeffs must be ordered the same in both of next 2 lines!
band_names = ([*coeffs.long_name] +
              [f'EOF{i}' for i in range(len(eofs))] +
              [f'EOF{i}_FOR_VIZ' for i in range(len(eofs_prepped))] +
              [*asynch.long_name]
             )
output = xr.concat([coeffs, eofs, eofs_prepped, asynch], dim='band')
assert len(band_names) == output.shape[0]
output = output.assign_coords({'band': band_names})

# write to disk
# NOTE: convert to xarray.Dataset so that band names write correctly to disk
output_filepath = os.path.join(phf.EXTERNAL_DATA_DIR,
                               'Terasaki_Hart_2024_LSP_coeffs_eofs_eofsforviz_asynch_FOR_GEE_DATA_VIEWER.tif',
                              )
output_ds = xr.Dataset()
for band_name in band_names:
    output_ds[band_name] = output.sel({'band': band_name})
output_ds.rio.to_raster(output_filepath)


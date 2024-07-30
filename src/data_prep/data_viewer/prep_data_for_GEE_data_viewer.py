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

# common SRS (unprojected lat,lon)
crs = 4326

# load coefficients
coeffs = rxr.open_rasterio(phf.COEFFS_FILE, masked=True)
assert coeffs.rio.crs == crs

# load asynchrony
neigh = 100
asynch = rxr.open_rasterio(phf.ASYNCH_FILES[neigh], masked=True)
assert asynch.rio.crs == crs

# load the raw EOFS
eofs = rxr.open_rasterio(phf.EOFS_FILE, masked=True)
eofs = eofs.rio.set_crs(crs)
assert eofs.rio.crs == crs

# load the scaled, ITCZ-folded EOFs, for data viz of phen significant MMRR results
eofs_prepped = rxr.open_rasterio(phf.EOFS_PREPPED_FILE)
eofs_prepped = eofs_prepped.rio.write_crs(8857)
eofs_prepped = eofs_prepped.rio.reproject_match(coeffs)

# stack them all and write to file
# NOTE: EOFs, asynch, then coeffs must be ordered the same in both of next 2 lines!
band_names = ([f'EOF{i}' for i in range(len(eofs))] +
              [f'EOF{i}_FOR_VIZ' for i in range(len(eofs_prepped))] +
              [*asynch.long_name] +
              [*coeffs.long_name]
             )
output = xr.concat([eofs, eofs_prepped, asynch, coeffs], dim='band')
assert len(band_names) == output.shape[0]
output = output.assign_coords({'band': band_names})

# write to disk
# NOTE: convert to xarray.Dataset so that band names write correctly to disk
output_filepath = os.path.join(phf.EXTERNAL_DATA_DIR,
                               'eofsraw_eofsfolded_asynch_coeffs_for_GEE_data_viewer.tif',
                              )
output_ds = xr.Dataset()
for band_name in band_names:
    output_ds[band_name] = output.sel({'band': band_name})
output_ds.rio.to_raster(output_filepath)



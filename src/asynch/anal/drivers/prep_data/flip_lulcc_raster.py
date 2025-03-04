import rioxarray as rxr
import os
import sys
sys.path.insert(1,
                '/global/home/users/drewhart/seasonality/seasonal_asynchrony/src/etc/')
import phen_helper_fns as phf
data_dir = phf.EXTERNAL_RF_DATA_DIR

flipped = rxr.open_rasterio(os.path.join(data_dir,
                                         'hansen_lulcc_pct_neigh_mean.tif'),
                            masked=True)
assert flipped.rio.crs.to_epsg() == 4326
unflipped = flipped.rio.reproject(4326)
unflipped.rio.to_raster(os.path.join(data_dir,
                                     'hansen_lulcc_pct_neigh_mean.tif'))

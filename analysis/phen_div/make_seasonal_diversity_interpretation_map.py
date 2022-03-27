import rioxarray as rxr
import numpy as np
import geopandas as gpd
import seaborn as sns
import colorsys
import math
import os

import helper_fns as hf

data_dir = '/home/deth/Desktop/CAL/research/projects/seasonality/results/maps'

# read in countries boundaries and EOFs
countries = gpd.read_file(os.path.join(data_dir, 'NewWorldFile_2020.shp'))
eofs = rxr.open_rasterio(os.path.join(data_dir, 'global_4_EOFs_coswts.tif'))
eofs_rotnorm = rxr.open_rasterio(os.path.join(data_dir,
                                    './global_4_EOFs_coswts_shemrot_normts.tif'))

# read in coeffs, for time-series fitting
coeffs = rxr.open_rasterio(os.path.join(data_dir, 'NIRv_global_coeffs.tif'))

# make design matrix
dm = hf.make_design_matrix()

def calc_time_series(coeffs, dm, normalize=False, lat_for_rotation=None):
    """
    calculate a time series from 5 coeffs
    """
    ts = np.sum(coeffs*dm, axis=1)
    if normalize:
        ts = (ts-np.min(ts))/(np.max(ts)-np.min(ts))
    if lat_for_rotation is not None:
        # TODO: DECIDE HOW TO ROTATE BASED ON LAT!
        pass
    return ts



# reproject to equal earth
countries = countries.to_crs(8857)
eofs=eofs.rio.set_crs('EPSG:4326')
eofs_rotnorm=eofs_rotnorm.rio.set_crs('EPSG:4326')
eofs_ee = eofs.rio.reproject("EPSG:8857")
eofs_rotnorm_ee = eofs_rotnorm.rio.reproject("EPSG:8857")

# plot both, in order
fig, ax_map = plt.subplots(1,1)
countries.plot(facecolor='none',
               edgecolor='black',
               linewidth=0.25,
               ax=ax_map,
              )
eofs_rotnorm_ee[:3,:,:].plot.imshow(ax=ax_map)

# create a circle-wheel, pick regularly spaced colors around it,
# sample cells' values within those colors, and calculate mean and std
# seasonality curves

# plot mean and std curves for each year (against winter solstice on x?)
# with some indicator of the color

# try clustering?

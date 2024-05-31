import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
from scipy.signal import find_peaks
from copy import deepcopy
import os, re

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf



# whether to run in debug mode
debug = False


#------------------------------------------------------------------------------
# fns
#------------------------------------------------------------------------------

def minmax_scale(vals, min_out=0, max_out=1):
    assert (min_out >=0) and (max_out <=1)
    scaled = (vals-np.min(vals))/(np.max(vals)-np.min(vals))
    scaled = scaled * (max_out - min_out)
    if min_out > 0:
        scaled = min_out + scaled
    return scaled


def rotate_time_series_to_min(ts,
                              nsteps=7,
                              return_cutidx=False,
                             ):
    '''
    shift time series so that it starts on its min value
    (or which of its min values is just before the greatest increase over the
    next n steps, if it has more than 1; e.g., zero-inflated observation data)
    '''
    if not isinstance(ts, np.ndarray):
        ts = np.array([*ts])
    minval = np.min(ts)
    #print('ts: ', ts)
    minidx = np.where(ts == minval)[0]
    #print('minidx: ', minidx)
    if len(minidx) == 1:
        #print('just 1')
        cutidx = minidx[0]
    else:
        #print('>1')
        tscat=  np.concatenate([ts]*2)
        minidx_nsteps_slope = (tscat[minidx+nsteps] - tscat[minidx])/nsteps
        #print('minidx_nsteps_slope: ', minidx_nsteps_slope)
        # NOTE: if there is more than one spot with the maximum observed
        #       n-step slope right after a min value then this will just default
        #       to the first one, which should be good enough
        minidx_maxdif = np.argmax(minidx_nsteps_slope)
        #print('minidx_maxdif: ', minidx_maxdif)
        cutidx = minidx[minidx_maxdif]
    #print('cutidx: ', cutidx)
    ts_rot = np.concatenate((ts[cutidx:], ts[:cutidx]))
    if return_cutidx:
        return ts_rot, cutidx
    else:
        return ts_rot


#------------------------------------------------------------------------------
# data load
#------------------------------------------------------------------------------

# load country boundaries
countries = gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                       'NewWorldFile_2020.shp')).to_crs(4326)

# load level-1 subnational jurisdictions (downloaded from:
#                                 https://gadm.org/download_country.html)
subnational = []
for f in [f for f in os.listdir(phf.BOUNDS_DIR) if re.search('^gadm.*json$', f)]:
    subnational.append(gpd.read_file(os.path.join(phf.BOUNDS_DIR,
                                                  f)).to_crs(4326))
subnational = pd.concat(subnational)

# load the coeffs
coeffs = rxr.open_rasterio(os.path.join(phf.EXTERNAL_DATA_DIR,
                                        '%s%s_coeffs.tif') % ('NIRv', ''))

# create the regression design matrix
# (for recreating the fitted LSP curves)
dm = phf.make_design_matrix()

# create array to store results
modality = np.ones(coeffs[0].values.shape) * np.nan

# get indices where coeffs are not null
I, J = np.where(pd.notnull(coeffs[0]))

# loop over pixels
for i, j in zip(I, J):
    # get the cell's fitted coefficients
    print(f"\n\tprocessing cell [{i}, {j}]...")
    coeffs_vec = coeffs[:, i, j].values
    # calc ts
    ts = phf.calc_time_series(coeffs_vec, dm)
    # minmax scale the ts
    mmts = minmax_scale(ts)
    # rotate to min value
    rmmts = rotate_time_series_to_min(mmts)
    # use simple neighbor-comparison & peak properties to calculate n peaks
    peaks, props = find_peaks(rmmts,
                             )
    if debug:
        fig, ax = plt.subplots(1, 1)
        ax.plot(rmmts, color='black')
        for peak in peaks:
            ax.axvline(peak, color='red', alpha=0.5)
        fig.show()
        input()
        plt.close('all')
    # calculate height ratio and store
    assert len(peaks) <= 2
    if len(peaks) == 2:
        heights = rmmts[peaks]
        diff = np.max(heights)-np.min(heights)
    else:
        diff = 1
    modality[i, j] = diff

# plot
modality_da = deepcopy(coeffs[0])
modality_da.values = modality
assert np.nanmin(modality) >= 0
assert np.nanmax(modality) == 1
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(1,1,1)
try:
    modality_da.plot.imshow(robust=True,
                            cmap='coolwarm_r',
                            ax=ax,
                            add_colorbar=True,
                            zorder=0,
                           )
except Exception:
    pass
subnational.plot(color='none',
                 edgecolor='black',
                 zorder=1,
                 ax=ax,
                 alpha=0.6,
                 linewidth=0.05,
                )
countries.plot(color='none',
               edgecolor='black',
               linewidth=0.2,
               zorder=2,
               ax=ax,
              )
ax.set_title('')
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_xticks(())
ax.set_yticks(())
fig.savefig('LSP_ann_sem_seasonality.png', dpi=600)

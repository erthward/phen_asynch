import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio as rio
from rasterio.plot import reshape_as_raster
import xarray as xr
import rioxarray as rxr
import os
import datetime
import daylight
import pytz
from shapely.geometry import Point
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from sklearn.manifold import MDS
from sklearn.preprocessing import normalize
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# TODO:
    # how to expand to whole raster?
    # could use clustering instead/as well, with scree plot to determine k


###################
# BEHAVIORAL PARAMS
###################

# size of random sample points to fit seasonality space
samp_size = 10_000

# data dir
data_dir = '/home/deth/Desktop/CAL/research/projects/seasonality/results/maps'

# set the seed
seed = 1
np.random.seed(seed)


###########
# FUNCTIONS
###########

def make_design_matrix():
    """
    Makes and returns the regression's design matrix, a 365 x 5 numpy array
    in which the columns contain, in order:
        - 1s (for the constant);
        - sin and cos of annual-harmonic days of the year
          (i.e. days are expressed in radians from 0 to 2pi);
        - sin and cos of the semiannual-harmonic days of the year
          (i.e. days are expressed in radians from 0 to 4pi).
    """
    # get 1 year of daily values, expressed in radians, 1 rotation/yr
    annual_radian_days = np.linspace(0, 2*np.pi, 366)[:365]
    # get 1 year of daily values, expressed in radians, 2 rotations/yr
    semiannual_radian_days = np.linspace(0, 4*np.pi, 366)[:365] % (2 * np.pi)
    # get the harmonic values of those
    sin1 = np.sin(annual_radian_days)
    cos1 = np.cos(annual_radian_days)
    sin2 = np.sin(semiannual_radian_days)
    cos2 = np.cos(semiannual_radian_days)
    # add a vector of 1s for the constant term, then recast as a 365 x 5 array,
    # to use as the covariate values in the regression
    design_mat = np.array([np.ones(sin1.shape), sin1, cos1, sin2, cos2]).T
    return design_mat


# functions to time-shift a 365-day time series to correct for latitude effect
# (for comparison across hemispheres)
def epoch(year, month, day, hour=0, minute=0, second=0, tz=pytz.UTC):
    return int(tz.localize(datetime.datetime(year, month, day, hour,
                                             minute, second)).timestamp())

def get_daylen_fractions(lat, return_all_days=False):

    # get sunclock object for this lat
    sun = daylight.Sunclock(lat, 0)
    # for each day of the year, subtract sunset from sunrise, divide by
    # daylength (in sec)
    year = 2022
    fullday = 60*60*24
    doys = []
    daylen_fracs = []

    for doy in range(1, 366):
        date = datetime.datetime.strptime(str(year) + "-" + str(doy), "%Y-%j")
        month = date.month
        day = date.day
        # NOTE: 2022 not a leap year
        sunrise = sun.sunrise(epoch(2022, month, day))
        sunset = sun.sunset(epoch(2022, month, day))
        diff = sunset-sunrise
        daylen_frac = diff/fullday
        doys.append(doy)
        daylen_fracs.append(daylen_frac)

    # find local maxima
    doys = np.array(doys)
    daylen_fracs = np.array(daylen_fracs)
    max_inds = argrelextrema(daylen_fracs, np.greater)[0]
    assert len(max_inds) <=2 or lat==0
    max_doys = doys[max_inds]
    # return doys of local maxima, and optionall full lists
    if return_all_days:
        return (max_doys, doys, daylen_fracs)
    else:
        return max_doys



#########################################
# READ IN COEFFS, GET RANDOM SAMPLE OF TS
#########################################

# read global NIRv coeffs
nirv = rxr.open_rasterio(os.path.join(data_dir, 'NIRv_global_coeffs.tif'))

# draw a random sample of points
X, Y = np.meshgrid(nirv.x, nirv.y)
samp = np.random.choice(range(X.size), size=samp_size, replace=False)
sampx = X.flatten()[samp]
sampy = Y.flatten()[samp]

# get coeffs at random points
samp_vals = np.stack([nirv.sel(x=sampx[i], y=sampy[i]) for i in range(len(sampx))])

# cast as df, remove NaNs
samp_df = pd.DataFrame(samp_vals)
samp_df.columns = nirv.long_name
samp_df['x'] = sampx
samp_df['y'] = sampy
samp_df = samp_df.dropna(how='any').reset_index().iloc[:, 1:]

# make design matrix
dm = make_design_matrix()

# get the time series for each sample
tss = np.ones((len(samp_df), 365)) * np.nan
for i, row in samp_df.iterrows():
    ts = np.sum(row.values[:5] * dm, axis=1)
    # rotate 1/2 year for soutern hemisphere
    # TODO: DOES THIS MAKE A DIFF? IT CHANGES COLORS, BUT SEEMS TO LEAVE
    # RELATIONSHIPS THE SAME!? NEED TO WORK THIS OUT NUMERICALLY
    if row['y']<0:
        ts = np.array([*ts[183:]] + [*ts[:183]])
    # normalize to [0,1]
    ts = normalize([ts]).flatten()
    assert ts.shape == (365,)
    tss[i, :] = ts





#############
# EMBED IN 3D
#############

# use metric multidimensional scaling to embed sample points' 365-length
# time series in 3-dimensional space (which I'll then map to RGB for mapping),
# using pairwise Euclidean distance between the 365-d points as dissim metric
mds = MDS(n_components=3, metric=True)
tss_3d = mds.fit_transform(tss)


#####################
# PLOT IN 3D, AND MAP
#####################

fig = plt.figure()
ax3d = fig.add_subplot(121, projection='3d')

# normalize the tss_3d columns, to then map onto RGB colors
tss_3d_norm = ((tss_3d - np.min(tss_3d, axis=0)) /
               (np.max(tss_3d, axis=0) - np.min(tss_3d, axis=0)))
# scatter on 3d axes and color by R (x), G (y), B (z) axes
ax3d.scatter(tss_3d[:, 0],
             tss_3d[:, 1],
             tss_3d[:, 2],
             c=tss_3d_norm,
             alpha=0.7)
ax3d.set_xlabel('red', fontdict={'fontsize':16})
ax3d.set_ylabel('green', fontdict={'fontsize':16})
ax3d.set_zlabel('blue', fontdict={'fontsize':16})


ax2 = fig.add_subplot(122)
countries = gpd.read_file(os.path.join(data_dir, 'NewWorldFile_2020.shp'))
countries = countries.to_crs(8857)
countries.plot(facecolor='none',
               edgecolor='black',
               linewidth=0.5,
               ax=ax2)
samp_gdf = gpd.GeoDataFrame(samp_df.drop(['x', 'y'], axis=1),
                            crs={'init': 'epsg:4326'},
                            geometry=[Point(xy) for xy in zip(samp_df.x,
                                                              samp_df.y)])
samp_gdf = samp_gdf.to_crs(8857)
ax2.scatter(samp_gdf.centroid.x,
            samp_gdf.centroid.y,
            c=tss_3d_norm, s=25,
            alpha=0.5)

# interpolate to grid
grid_tss_3d = griddata(samp_df.loc[:, ['x', 'y']].values, tss_3d_norm,
                       (X, Y), method='cubic')

# mask out NAs
for i in range(3):
    grid_tss_3d[:,:,i][np.isnan(nirv[0])] = np.nan

# write to file
with rio.Env():
    profile = rio.open(os.path.join(data_dir, 'NIRv_global_coeffs.tif')).profile
    profile.update(count=3, compress='lzw')
    with rio.open('seasonality_MDS_result.tif', 'w', **profile) as dst:
        # NOTE: profile from input file specifies float32, originating from GEE
        dst.write(reshape_as_raster(grid_tss_3d).astype(rio.float32))

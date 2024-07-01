import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
from shapely.geometry import Point
from shapely.geometry import Polygon as shapelyPolygon
from matplotlib.collections import PatchCollection
import numpy as np
import xarray as xr
import rasterio as rio
import rioxarray as rxr
from datetime import datetime
import statsmodels.api as sm
import zipfile36 as zipfile
import seaborn as sns
import re, os, sys

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                                                            'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf



# plotting params
axlabel_fontdict={'fontsize': 14}
ticklabel_fontsize=10

# main behavioral params
rescale=True
delete_after_finished = True
plot_time_series = True
max_neigh_cell_dist = 2 # at our 0.05deg res, this is up to ~11km away...
seed = 1
np.random.seed(seed)

# choose dataset
rs_var = sys.argv[1]
assert rs_var in ['SIF', 'NIRv']

# choose masking mod
masking_mode = 'default'
#masking_mode = 'strict'
masking_suffix = '_STRICT' * (masking_mode == 'strict')


rs_var_units = {'SIF': '$mW\ m^{-2}\ sr^{-1}\ nm^{-1}$',
                       'NIRv': '$nmol\ m^{-2}\ s^{-1}$'
               }


filter_start_date = None
filter_end_date = None


# data directories
rs_datadir = phf.EXTERNAL_DATA_DIR
flux_datadir = phf.EXTERNAL_FLUX_DATA_DIR
other_datadir = phf.DATA_DIR

# indicate the RS-based coefficients TIFF to evaluate
rs_coeffs_tif = os.path.join(rs_datadir,
                             '%s_coeffs%s.tif' % (rs_var, masking_suffix))

# variables that differ between main FLUXNET sites and CH4 sites
# (indexed by the filename patterns that distinguish between the two types)
response_vars = {'CH4': 'GPP_DT',
                 'SUBSET': 'GPP_DT_VUT_50',
                }
qc_cols = {'CH4': '',
           'SUBSET': 'NEE_VUT_50_QC',
          }

# load table with all fluxent sites' coordinates and info
sites = pd.read_csv('./FLUXNET_sites.csv')

# load the RS-based seasonality coeffs file, to be evaluated
def load_rs_coeffs(rs_coeffs_tif):
    rs_coeffs = rxr.open_rasterio(rs_coeffs_tif)
    # NOTE: already in WGS84 latlon, and I'll just be extracting using
    #       latlon, so no need to transform
    return rs_coeffs


# get site info corresponding to a filename
def get_site_info(fn):
    # get rid of path string
    fn = os.path.split(fn)[-1]
    # patt to get site id from filename
    patt = '(?<=^FLX_)\w{2}-\w{3}(?=_FLUXNET)'
    # get the site ID
    id = re.search(patt, fn).group()
    # look up this ID's row in sites table
    row = sites[sites.SITE_ID == id]
    assert len(row) == 1
    row_info = row.iloc[0]
    name = row_info['SITE_NAME']
    lon = row_info['LOCATION_LONG']
    lat = row_info['LOCATION_LAT']
    igbp = row_info['IGBP']
    mat = row_info['MAT']
    map = row_info['MAP']
    return id, name, lon, lat, igbp, mat, map


# fit the harmonic regression to each, then
# get the detrended, predicted values for each site
def fit_harmonic_regression(df, response):
    # gather predictors into an array, including:
    # day of total time series (to detrend)
    t = phf.minmax_rescale_array(df.index)
    # sines and cosines of annual and semiannual components
    sin_ann = [np.sin(n) for n in df['ann']]
    cos_ann = [np.cos(n) for n in df['ann']]
    sin_sem = [np.sin(n) for n in df['sem']]
    cos_sem = [np.cos(n) for n in df['sem']]
    X = np.stack((t, sin_ann, cos_ann, sin_sem, cos_sem)).T
    # and grab the y
    y = np.atleast_2d(np.array(df[response])).T
    # build and fit the regression
    reg = sm.OLS(y, X).fit()
    return reg


# little algo to create a nested array of ring dists
def get_nested_array_of_ring_dists(max_dist):
    arr = np.ones([max_dist*2+1]*2)*max_dist
    for m in range(1,max_dist):
        arr[m:-m, m:-m] = max_dist-m
    arr[max_dist, max_dist] = 0
    return arr


def get_all_nearby_cells(x, y, cell_res, max_neigh_cell_dist):

    # get a series of values stepping max_neigh_cell_dist res steps
    # in each direction from the focal x,y
    nearby_x = np.linspace(x-(max_neigh_cell_dist*cell_res),
                            x+(max_neigh_cell_dist*cell_res),
                            2*max_neigh_cell_dist+1)
    nearby_y = np.linspace(y-(max_neigh_cell_dist*cell_res),
                            y+(max_neigh_cell_dist*cell_res),
                            2*max_neigh_cell_dist+1)

    # get all coord pairs for those nearby locations
    nearby_X, nearby_Y = np.meshgrid(nearby_x, nearby_y)

    # create a square array of neigh_cell ring distances
    cell_dists = get_nested_array_of_ring_dists(max_neigh_cell_dist)

    # melt everything down and return it
    nearby_X = nearby_X.ravel()
    nearby_Y = nearby_Y.ravel()
    cell_dists = cell_dists.ravel()
    nearbys = [*zip(nearby_X, nearby_Y)]

    # dict to return results keyed by increasing cell dist
    nearby_cells_by_dist = {}
    for curr_dist in range(1, max_neigh_cell_dist+1):
        nearbys_curr_dist = [nb for nb, cd in zip(nearbys, cell_dists) if cd ==
                             curr_dist]
        # shuffle them
        np.random.shuffle(nearbys_curr_dist)
        nearby_cells_by_dist[curr_dist] = nearbys_curr_dist

    # convert that to a list to popped through
    nearby_cells_to_check = []
    for cell_dist in sorted([*nearby_cells_by_dist.keys()]):
        cells_to_check_this_dist = [*zip([cell_dist]*len(
                                         nearby_cells_by_dist[cell_dist]),
                                         nearby_cells_by_dist[cell_dist])]
        nearby_cells_to_check.extend(cells_to_check_this_dist)

    # invert list, so that list.pop() draws from the closest cells outward
    nearby_cells_to_check = nearby_cells_to_check[::-1]

    return nearby_cells_to_check



def predict_rs_detrended_vals(coeffs_rast, x, y, design_mat,
                              max_neigh_cell_dist, rescale=True):
    """
    Calculates the predicted time series at pixel at x,y in a rioxarray raster,
    using the coefficients for the constant and the
    sin and cosine terms of the annual and semiannual
    harmonic components. Returns the time series as a numpy array.
    """
    # get the nearest pixel's coefficients from the rioxarray dataset
    coeffs = coeffs_rast.sel(x=x, y=y, method="nearest").values
    cell_dist = 0

    # if coeffs are nans then try at cells out to max_neigh_cell_dist
    if np.all(pd.isnull(coeffs)):
        cell_res = float(coeffs_rast.spatial_ref.GeoTransform.split()[1])
        nearby_cells = get_all_nearby_cells(x, y, cell_res,
                                              max_neigh_cell_dist)
        # step through increasing neigh cell rings, looking for one with non-nan
        # coeff vals to use
        while np.all(pd.isnull(coeffs)) and len(nearby_cells) > 0:
            next_cell_dist, next_cell_coords = nearby_cells.pop()
            cell_dist = next_cell_dist
            x,y = next_cell_coords
            coeffs = coeffs_rast.sel(x=x, y=y, method="nearest").values

    # if failed to find coeffs, set cell_dist to np.nan
    if np.all(pd.isnull(coeffs)):
        cell_dist = np.nan


    # multiply the pixel's set of coefficients by the design mat, then sum
    # all the regression terms
    # NOTE: the coeffs are a numpy array of shape (5,);
    #       the design matrix is a numpy array of shape (365, 5);
    #       the multiplication operates row by row, and elementwise on each
    #       row, returning a new (365, 5) array that, when summed across the
    #       columns (i.e. the regression's linear terms) gives a (365,) vector
    #       of fitted daily values for the pixel
    pred = np.sum(coeffs * design_mat, axis=1)

    # rescale it
    if rescale:
        pred = phf.minmax_rescale_array(pred)

    # pair it with a date range object, then make into a df
    dates = pd.date_range(start='1/1/2021', end='12/31/2021')
    pred_df = pd.DataFrame({'date': dates, 'rs_pred': pred})
    return pred_df, cell_dist


def calc_euc_dist(a1, a2):
    """
    Calculates the Euclidean distance between two 1d, length-n numpy arrays.

    Returns the distance as a float.
    """
    dist = np.sqrt(np.sum((a1 - a2)**2))
    return dist


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


# predict the detrended values
def predict_fluxnet_detrended_vals(df, mod, rescale=True):
    # sines and cosines of annual and semiannual components
    sin_ann = [np.sin(n) for n in df['ann']]
    cos_ann = [np.cos(n) for n in df['ann']]
    sin_sem = [np.sin(n) for n in df['sem']]
    cos_sem = [np.cos(n) for n in df['sem']]
    X = np.stack((np.ones(len(sin_ann)), sin_ann, cos_ann, sin_sem, cos_sem)).T
    predicted = mod.predict(X)
    # rescale it
    if rescale:
        predicted = phf.minmax_rescale_array(predicted)
    df['pred'] = predicted


# process the site's data and return the processed df
def process_site_data(zip_filename,
                      rescale=True,
                      filter_start_date=None,
                      filter_end_date=None,
                      delete_after_finished=False):

    # determine if file is main FLUXNET or CH4 site
    sitetype_patt = re.search('(SUBSET)|(CH4)', zip_filename).group()

    # get response variable, based on sitetype_patt
    response_var = response_vars[sitetype_patt]

    # unzip the archive
    with zipfile.ZipFile(zip_filename, 'r') as z:
        # get the right file to extract
        filelist = z.filelist
        csv_patt = re.compile('%s_DD' % sitetype_patt)
        potential_csv_file = [f for f in filelist if re.search(csv_patt,
                                                               f.filename)]
        assert len(potential_csv_file) == 1
        csv_file = potential_csv_file[0]
        # extract that file to cwd
        csv_filename = z.extract(csv_file)

    # read the CSV
    df = pd.read_csv(csv_filename)

    # drop data with mean daily QC val <0.7
    # (per Reviewer 3 recommendation)
    if sitetype_patt == 'SUBSET':
        df = df[df.loc[:, qc_cols[sitetype_patt]] >= 0.7]

    # make datestamp columns in datetime objects
    df['date'] = pd.to_datetime(df.loc[:, 'TIMESTAMP'], format='%Y%m%d')

    # subset GPP vars
    gpp = df.loc[:,[*df.columns[[col.startswith('GPP') for col in
                                         df.columns]]]+['date']]
    gpp.set_index('date', inplace=True)

    # drop missing data
    gpp = gpp[gpp[response_var] != -9999]

    # filter out dates, if needed
    if filter_start_date is not None and filter_end_date is not None:
        # for now, just drop years after 2006 in the forest dataset,
        # because of the big gap caused by the missing 2007-ish data
        # TODO: COME UP WITH BETTER SOLUTION
        gpp = gpp.loc[filter_start_date:filter_end_date]

    # add numeric day of year columns,
    # and add annual and semi-annual circular time compononet columns (radians)
    gpp['doy'] = [idx.timetuple().tm_yday for idx in gpp.index]
    gpp['ann'] = 2*np.pi*np.array([*gpp['doy']])/365
    gpp['sem'] = 2*np.pi*np.array([*gpp['doy']])/(365/2)

    # save length of GPP time series being used for the regression (in years)
    gpp_ts_len = (gpp.index[-1] - gpp.index[0]).days/365.25

    # fit the regression, then add the predicted, detrended vals col
    reg = fit_harmonic_regression(gpp, response_var)
    predict_fluxnet_detrended_vals(gpp, reg, rescale=rescale)

    # get the overall model P-value
    pval = reg.f_pvalue

    # delete the file, if requested
    if delete_after_finished:
        os.remove(csv_filename)

    return gpp, gpp_ts_len, pval


# get and compare flux-based and RS-based predicted GPP values for a site
def compare_rs_flux_predicted_vals(zip_filename, coeffs_rast, design_mat,
                                   max_neigh_cell_dist,
                                   rescale=True,
                                   filter_start_date=None,
                                   filter_end_date=None,
                                   delete_after_finished=False,
                                   plot_time_series=False):
    # get site info
    id, name, lon, lat, igbp, mat, map = get_site_info(zip_filename)

    # get rs-predicted seasonality
    rs_pred, cell_dist = predict_rs_detrended_vals(coeffs_rast, lon, lat,
                                                   design_mat,
                                                   max_neigh_cell_dist,
                                                   rescale=rescale)

    # return NaN, if this flux tower is not located within striking distance
    # of a valid seasonality-fitted pixel
    if np.all(pd.isnull(rs_pred['rs_pred'])):
        dist = np.nan
        r2 = np.nan
        flux_pred_df = rs_pred = None
        gpp_ts_len = np.nan
        gpp_pval = np.nan

    else:

        # get df with fluxnet-predicted seasonality,
        # and length of GPP time series being used to fit the regression
        flux_pred_df, gpp_ts_len, gpp_pval = process_site_data(zip_filename,
                                    rescale=rescale,
                                    filter_start_date=filter_start_date,
                                    filter_end_date=filter_end_date,
                                    delete_after_finished=delete_after_finished,
                                                              )
        # 'rotate' that data to fill a 'standard' year (01/01 to 12/31, 2021),
        # or at least as much of that year as we have days of the year to fill
        dates = pd.date_range(start='1/1/2021', end='12/31/2021')
        date_preds = []
        for date in dates:
            doy = date.day_of_year
            subdf = flux_pred_df[flux_pred_df['doy'] == doy]
            if len(subdf) > 0:
                date_pred = subdf.iloc[0,:]['pred']
            else:
                date_pred = np.nan
            date_preds.append(date_pred)
        flux_pred = pd.DataFrame({'date': dates, 'flux_pred': date_preds})

        # merge both to get matched predictions for all common dates
        merged = pd.merge(rs_pred, flux_pred, how='inner', on='date')

        # drop NaNs
        merged_nonans = merged.dropna()

        # calculate the Euclidean distance between both time series
        dist = calc_euc_dist(merged_nonans['rs_pred'],
                             merged_nonans['flux_pred'],
                            )

        # calculate the R^2 between both time series
        r2 = np.corrcoef(merged_nonans['rs_pred'],
                          merged_nonans['flux_pred'],
                        )[0,1]**2

        # plot and save, if requested
        if plot_time_series:
            fig, ax = plt.subplots(1)
            ax.plot(merged['rs_pred'], '-k', label=rs_var)
            ax.plot(merged['flux_pred'], '-r', label='FLUX')
            ax.set_xlabel("day of year", fontdict={'fontsize':9})
            ax.set_ylabel(("rescaled metric of seasonality\n"
                           "'%s'=LSP (%s); "
                           #"'FLUX'=GPP ($\mu mol\ CO_2\ m^{-2}\ s{-1}$"
                           "'FLUX'=GPP ($g\ C\ m^{-2}\ d{-1}$"
                          ")") % (rs_var, rs_var_units[rs_var]),
                          fontdict={'fontsize':9})
            ax.legend()
            fig.suptitle('%s: %s: DIST=%0.3f; $R^2$=%0.3f' % (id, name,
                                                              dist, r2))
            fig.savefig(os.path.join(phf.FIGS_DIR, './flux_plots/%s_%s.png' % (id,
                                                                          rs_var)))
            plt.close('all')


    # gather summary data into a dict, for ingestion into a summary df
    result_dict = {'id': id,
                   'name': name,
                   'lon': lon,
                   'lat': lat,
                   'igbp': igbp,
                   'mat': mat,
                   'map': map,
                   'dist': dist,
                   'r2': r2,
                   'cell_dist': cell_dist,
                   'gpp_ts_len': gpp_ts_len,
                   'gpp_pval': gpp_pval,
                   'notes': np.nan,
                  }

    return result_dict, rs_pred, flux_pred_df


# load the rs-derived coefficients (fitted on GEE)
coeffs_rast = load_rs_coeffs(rs_coeffs_tif)

# make design matrix used to estimate rs-based coeffs' fitted seasonality
design_mat = make_design_matrix()


# loop over all zipfiles, run evaluation, save evaluation metrics
results = {'id': [],
           'name': [],
           'lon': [],
           'lat': [],
           'igbp': [],
           'mat': [],
           'map': [],
           'dist': [],
           'r2': [],
           'cell_dist': [],
           'gpp_ts_len': [],
           'gpp_pval': [],
           'notes': [],
          }

zip_filenames = [f for f in os.listdir(flux_datadir) if
                 os.path.splitext(f)[-1] == '.zip']
zip_filenames = [os.path.join(flux_datadir, fn) for fn in zip_filenames]

# NOTE: decided to skip CH4 sites altogether because a number of them I've
#       looked at have no clear QC flags included for the GEE data
zip_filenames = [f for f in zip_filenames if re.search('FLUXNET-CH4_', f) is None]

for zip_filename in zip_filenames:
    try:
        print('\n\nNow processing: %s\n\n' % zip_filename)
        result, rs_pred, flux_pred_df = compare_rs_flux_predicted_vals(
                                   zip_filename, coeffs_rast, design_mat,
                                   max_neigh_cell_dist,
                                   rescale=rescale,
                                   filter_start_date=filter_start_date,
                                   filter_end_date=filter_end_date,
                                   delete_after_finished=delete_after_finished,
                                   plot_time_series=plot_time_series)
    except Exception as e:
        id, name, lon, lat, igbp, mat, map = get_site_info(zip_filename)
        print('\n\nSite %s: %s failed with error:\n\n\t%s\n\n' % (id, name, e))
        result = {'id': id,
                  'name': name,
                  'lon': lon,
                  'lat': lat,
                  'igbp': igbp,
                  'mat': mat,
                  'map': map,
                  'dist': np.nan,
                  'r2': np.nan,
                  'cell_dist': np.nan,
                  'gpp_ts_len': np.nan,
                  'gpp_pval': np.nan,
                  'notes': e
                 }

    for k,v in result.items():
        results[k].append(v)

    print('-'*80)

# convert results to a dataframe
results_df = pd.DataFrame(results)

# drop datasets with duplicates at the same site
# (1 from main FLUXNET, 1 from CH4)
# because of strong non-independence),
# keeping whichever one has the longer GPP time series available
results_df = results_df.sort_values(['id', 'gpp_ts_len'],
            ascending=[True, True]).drop_duplicates(subset='id',
                                                    keep='last')
# and drop sites labeled as water bodies, urban, and cropland
results_df = results_df[~results_df['igbp'].isin(['WAT', 'URB', 'CRO'])]

# save all results
results_df.to_csv('FLUXNET_evaluation_results_%s.csv' % rs_var, index=False)


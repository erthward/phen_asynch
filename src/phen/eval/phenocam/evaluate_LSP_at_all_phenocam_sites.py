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
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf



# plotting params
axlabel_fontdict={'fontsize': 14}
ticklabel_fontsize=10

# main behavioral params
rescale=True
plot_time_series = True
max_neigh_cell_dist = 2 # at our 0.05deg res, this is up to ~11km away...
seed = 1
np.random.seed(seed)

# choose dataset
try:
    rs_var = sys.argv[1]
except Exception:
    rs_var = input('Which RS var?...')
assert rs_var in ['SIF', 'NIRv']

# choose masking mod
masking_mode = 'default' # or 'strict'
masking_suffix = '_STRICT' * (masking_mode == 'strict')


rs_var_units = {'SIF': '$mW\ m^{-2}\ sr^{-1}\ nm^{-1}$',
                       'NIRv': '$nmol\ m^{-2}\ s^{-1}$'
               }


# data directories
rs_datadir = phf.EXTERNAL_DATA_DIR
phenocam_datadir = phf.EXTERNAL_PHENOCAM_DATA_DIR

# indicate the RS-based coefficients TIFF to evaluate
rs_coeffs_tif = os.path.join(rs_datadir,
                             '%s_coeffs%s.tif' % (rs_var, masking_suffix))

# response var to model in PhenoCam datasets
# NOTE: 75th percentile of GCC best balances the desire to model the
#       mean/center value with Richardson et al. 2018 comment
#       that they generally find 90th percentile best tracks the mode without
#       sensitivity to lighting conditions but occasionally tracks extreme
#       values too closely
#response_var='gcc_75'
response_var = 'ndvi_75'

# load table with all PhenoCam ROIs' coordinates and info
rois_df = pd.read_csv('ROIs.csv')

# load table with site-level metadata
sites_df = pd.read_csv('site_metadata.csv')

# load the RS-based seasonality coeffs file, to be evaluated
def load_rs_coeffs(rs_coeffs_tif):
    rs_coeffs = rxr.open_rasterio(rs_coeffs_tif)
    # NOTE: already in WGS84 latlon, and I'll just be extracting using
    #       latlon, so no need to transform
    return rs_coeffs


# get ROI info corresponding to a filename
def get_roi_info(filepaths, rois_df):
    # use first filepath for site info
    filepath = filepaths[0]
    # get rid of path string
    filename = os.path.split(filepath)[-1]
    # get the site ID
    site = re.search("[a-z,A-Z,0-9,\.,\-]+(?=_)", os.path.split(filename)[-1]).group()
    vegs = []
    files_3day = []
    for fp in filepaths:
        fn = os.path.split(fp)[-1]
        veg = re.search("(?<=_)[A-Z]{2}(?=_)", os.path.split(fn)[-1]).group()
        vegs.append(veg)
        roi_id = re.search("(?<=_)\d{1,4}(?=_" + f"{response_var.split('_')[0]}_3day\.csv)",
                           os.path.split(fn)[-1]).group()
        roi_name = f"{site}_{veg}_{roi_id.zfill(4)}"
        # look up this ID's row in sites table
        row = rois_df[rois_df['roi_name'] == roi_name]
        assert len(row) == 1
        row_info = row.iloc[0]
        lon = row_info['lon']
        lat = row_info['lat']
        file_3day = row_info['three_day_summary']
        files_3day.append(file_3day)
    return site, lon, lat, vegs, files_3day


def read_phenocam_ndvi_file(filepath):
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            header_line_nums = [i for i, l in enumerate(lines) if re.search('^date,year,doy', l)]
            assert len(header_line_nums) == 1
            header_line_num = header_line_nums[0]
        df = pd.read_csv(filepath, skiprows=header_line_num)
    except Exception:
        try:
            df = pd.read_csv(filepath)
        except Exception:
            df = None
    return df


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



def predict_detrended_lsp_vals(coeffs_rast,
                               x,
                               y,
                               design_mat,
                               max_neigh_cell_dist,
                               rescale=True,
                              ):
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
def predict_phenocam_detrended_vals(df, mod):
    # sines and cosines of annual and semiannual components
    sin_ann = [np.sin(n) for n in df['ann']]
    cos_ann = [np.cos(n) for n in df['ann']]
    sin_sem = [np.sin(n) for n in df['sem']]
    cos_sem = [np.cos(n) for n in df['sem']]
    X = np.stack((np.ones(len(sin_ann)), sin_ann, cos_ann, sin_sem, cos_sem)).T
    predicted = mod.predict(X)
    df['pred_raw'] = predicted
    df['pred'] = phf.minmax_rescale_array(predicted)


# process the ROI's data and return the processed df
def process_roi_data(filepaths_list,
                     response_var=response_var,
                     min_req_len_yrs=1,
                    ):
    # get file's site ID
    site = re.search("[a-z,A-Z,0-9,\.,\-]+(?=_)", os.path.split(filepaths_list[0])[-1]).group()

    # read the data
    dfs = []
    for filepath in filepaths_list:
        if 'ndvi' in response_var:
            dfs.append(read_phenocam_ndvi_file(filepath))
        else:
            dfs.append(pd.read_csv(filepath))
    df = pd.concat(dfs)

    # checks on structure
    assert response_var in df.columns
    assert 'doy' in df.columns

    # make datestamp columns in datetime objects
    df['date'] = pd.to_datetime(df.loc[:, 'date'], format='%Y-%m-%d')
    df.set_index('date', inplace=True)

    # subset to key columns
    outlier_var = 'outlierflag_'+response_var
    df = df.loc[:, ['doy', response_var, outlier_var]]

    # drop outliers
    non_outliers = (pd.isnull(df[outlier_var])) | (df[outlier_var]==0)
    df = df[non_outliers]

    # drop missing data
    df = df[pd.notnull(df[response_var])]

    # make sure not length zero, otherwise return NaNs
    if len(df) == 0:
        raise ValueError((f"time series for {site} doesn't satisfy "
                          "minimum length!"))

    # sort by date
    df = df.sort_index()

    # get length of GPP time series being used for the regression (in years)
    ts_len = (df.index[-1] - df.index[0]).days/365.25

    # check if we meet the minimum required time series length
    if not ts_len >= min_req_len_yrs:
        raise ValueError((f"time series for {site} doesn't satisfy "
                          "minimum length!"))

    # add annual and semi-annual circular time compononet columns (radians)
    df['ann'] = 2*np.pi*np.array([*df['doy']])/365
    df['sem'] = 2*np.pi*np.array([*df['doy']])/(365/2)

    # average different ROIs' values by date
    df = df.groupby(df.index.date).mean()

    # convert doy back to integer
    df['doy'] = [int(doy) for doy in df['doy'].values]

    # fit the regression, then add the predicted, detrended vals col
    reg = fit_harmonic_regression(df, response_var)
    predict_phenocam_detrended_vals(df, reg)

    # get the model's R^2 and overall P-value
    r2 = reg.rsquared
    pval = reg.f_pvalue

    return df, ts_len, r2, pval


# get and compare LSP and PhenoCam values
def compare_predicted_vals(filepaths_list,
                           rois_df,
                           coeffs_rast,
                           design_mat,
                           max_neigh_cell_dist,
                           rescale=True,
                           plot_time_series=False,
                          ):
    # get site info
    site, lon, lat, veg, files_3day = get_roi_info(filepaths_list,
                                                   rois_df,
                                                  )

    # get rs-predicted seasonality
    rs_pred, cell_dist = predict_detrended_lsp_vals(coeffs_rast,
                                                    lon,
                                                    lat,
                                                    design_mat,
                                                    max_neigh_cell_dist,
                                                    rescale=rescale,
                                                   )

    # return NaN, if this PhenoCam site is not located within striking distance
    # of a valid seasonality-fitted pixel
    if np.all(pd.isnull(rs_pred['rs_pred'])):
        dist = np.nan
        cam_reg_r2 = np.nan
        r2 = np.nan
        cam_pred_df = rs_pred = None
        pval = np.nan
        ts_len = np.nan

    else:

        # get df with raw and fitted PhenoCam GCC values, PhenoCam
        # time series length, and PhenoCam regression P-value
        cam_pred_df, ts_len, cam_reg_r2, pval = process_roi_data(filepaths_list)

        # 'rotate' that data to fill a 'standard' year (01/01 to 12/31, 2021),
        # or at least as much of that year as we have days of the year to fill
        dates = pd.date_range(start='1/1/2021', end='12/31/2021')
        date_preds = []
        for date in dates:
            doy = date.day_of_year
            subdf = cam_pred_df[cam_pred_df['doy'] == doy]
            if len(subdf) > 0:
                date_pred = subdf.iloc[0,:]['pred']
            else:
                date_pred = np.nan
            date_preds.append(date_pred)
        cam_pred = pd.DataFrame({'date': dates,
                                 'cam_pred': date_preds,
                                })

        # merge both to get matched predictions for all common dates
        merged = pd.merge(rs_pred, cam_pred, how='inner', on='date')

        # drop NaNs
        merged_nonans = merged.dropna()

        # calculate the Euclidean distance between both time series
        dist = calc_euc_dist(merged_nonans['rs_pred'],
                             merged_nonans['cam_pred'],
                            )

        # calculate the R^2 between both time series
        r2 = np.corrcoef(merged_nonans['rs_pred'],
                         merged_nonans['cam_pred'],
                        )[0,1]**2

        # plot and save, if requested
        if plot_time_series:
            fig, axs = plt.subplots(2,1)
            ax = axs[0]
            ax.plot(merged['rs_pred'], '.k', label=rs_var)
            ax.plot(merged['cam_pred'].values, '.r', label='PhenoCam')
            ax.set_xlabel("day of year", fontdict={'fontsize':9})
            ax.set_ylabel(("rescaled metric of seasonality\n"
                           "'%s'=LSP (%s); "
                           "'PhenoCam'=GCC %s") % (rs_var,
                                                   rs_var_units[rs_var],
                                                   response_var.split('_')[-1],
                                                  ),
                          fontdict={'fontsize':9})
            ax.legend()
            ax = axs[1]
            ax.plot(cam_pred_df.index, cam_pred_df[response_var], '.k', zorder=1)
            for yr in np.unique([d.year for d in cam_pred_df.index]):
                ax.axvline(pd.Timestamp(datetime(yr, 1, 1, 0, 0, 0)),
                           -1,
                           1,
                           color='gray',
                           linestyle=':',
                           alpha=0.5,
                           zorder=0,
                          )
            ax.set_xlabel('date')
            ax.set_ylabel(response_var)
            ax.tick_params(labelsize=6, rotation=35)
            fig.suptitle('%s: PhenoCam reg $R^2$=%0.3f; PhenoCam-%s$R^2$=%0.3f' % (site,
                                                                                   cam_reg_r2,
                                                                                   rs_var,
                                                                                   r2,
                                                                                  ))
            fig.savefig(os.path.join(phf.FIGS_DIR,
                        'phenocam_ts_plots/%s_%s.png' % (site, rs_var)))
            plt.close('all')


    # gather summary data into a dict, for ingestion into a summary df
    result_dict = {'site': site,
                   'lon': lon,
                   'lat': lat,
                   'veg': ';'.join(veg),
                   'dist': dist,
                   'cam_reg_r2': cam_reg_r2,
                   'r2': r2,
                   'cell_dist': cell_dist,
                   'cam_ts_len': ts_len,
                   'pval': pval,
                   'files': ';'.join(files_3day),
                   'notes': np.nan,
                  }

    return result_dict, rs_pred, cam_pred_df


# load the rs-derived coefficients (fitted on GEE)
coeffs_rast = load_rs_coeffs(rs_coeffs_tif)

# make design matrix used to estimate rs-based coeffs' fitted seasonality
design_mat = make_design_matrix()

# loop over all zipfiles, run evaluation, save evaluation metrics
results = {'site': [],
           'lon': [],
           'lat': [],
           'veg': [],
           'dist': [],
           'cam_reg_r2': [],
           'r2': [],
           'cell_dist': [],
           'cam_ts_len': [],
           'pval': [],
           'files': [],
           'notes': [],
          }

filenames = [f for f in os.listdir(phenocam_datadir) if
                 os.path.splitext(f)[-1] == '.csv']
# drop anything that isn't a PhenoCam dataset (e.g., results tables)
datafile_patt = ("[a-z,A-Z,0-9,\.,\-]_[A-Z]{2}_\d{1,4}_"
                 f"{response_var.split('_')[0]}_3day.csv")
filenames = [f for f in filenames if re.search(datafile_patt, f) is not None]

filepaths = [os.path.join(phenocam_datadir, fn) for fn in filenames]

sites_dict = {}
for fp in filepaths:
    site = os.path.split(fp)[1].split('_')[0]
    if site in sites_dict:
        sites_dict[site].append(fp)
    else:
        sites_dict[site] = [fp]

for site, filepaths_list in sites_dict.items():
    try:
        print('\n\nNow processing %s comparison for site %s\n\n' % (rs_var, site))
        (result,
         rs_pred,
         cam_pred_df) = compare_predicted_vals(filepaths_list,
                                               rois_df,
                                               coeffs_rast,
                                               design_mat,
                                               max_neigh_cell_dist,
                                               rescale=rescale,
                                               plot_time_series=plot_time_series,
                                              )
    except ValueError as e:
        site, lon, lat, veg, files_3day = get_roi_info(filepaths_list,
                                                      rois_df,
                                                     )
        print('\n\nSite %s failed with error:\n\n\t%s\n\n' % (site, e))
        result = {'site': site,
                  'lon': lon,
                  'lat': lat,
                  'veg': ';'.join(veg),
                  'dist': np.nan,
                  'cam_reg_r2': np.nan,
                  'r2': np.nan,
                  'cell_dist': np.nan,
                  'cam_ts_len': np.nan,
                  'pval': np.nan,
                  'files': ';'.join(files_3day),
                  'notes': e
                 }

    for k,v in result.items():
        results[k].append(v)

    print('-'*80)

# convert results to a dataframe
results_df = pd.DataFrame(results)

# merge in site metadata
sites_subdf = sites_df.loc[:, ['sitename',
                               'long_name',
                               'elevation',
                               'date_start',
                               'date_end',
                               'MAP_worldclim',
                               'MAT_worldclim',
                               'primary_veg_type',
                               'secondary_veg_type',
                               'landcover_igbp',
                               'site_acknowledgements',
                              ]]
results_df_aug = pd.merge(results_df,
                          sites_subdf,
                          left_on='site',
                          right_on='sitename',
                          how='inner',
                         )
assert len(results_df_aug) == len(results_df)

# drop sites with invalid IGBP land cover and with ag
has_ag = np.array(['AG' in veg for veg in results_df_aug['veg'].values])
has_ag = (has_ag |
          ([(pd.notnull(v) and
             v == 'AG') for v in results_df_aug['primary_veg_type'].values]))
has_ag = (has_ag |
          ([(pd.notnull(v) and
             v == 'AG') for v in results_df_aug['secondary_veg_type'].values]))
# NOTE: per PhenoCam site, and matching docs elsewhere, these values code
#       for all forest, shrubland, savanna, grassland, and permanent wetland
# NOTE: just assume missing IGBP values are valid land cover, since this isn't
#       a mission critical decision anyhow
valid_lc_classes = [1,2,3,4,5,6,7,8,9,10,11]
valid_lc = (~has_ag & (results_df_aug['landcover_igbp'].isin(valid_lc_classes) |
                       pd.isnull(results_df_aug['landcover_igbp'])))
results_filt = results_df_aug[valid_lc]

# save all results
results_filt.to_csv('PhenoCam_evaluation_results_%s.csv' % rs_var, index=False)



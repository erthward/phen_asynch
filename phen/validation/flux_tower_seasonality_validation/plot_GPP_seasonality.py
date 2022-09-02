import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score


# TODO:
# use filename to grab lat and lon from fluxnet_sites.csv,
# then use that to grab sat-fitted seasonality and plot that as well!


# overarching params
response_fluxnet = 'GPP_DT_VUT_50'
response_ch4 = 'GPP_DT'



# get full-year, contiguous slices of data
def keep_only_full_years(df):
    yrs = np.unique([str(idx)[:4] for idx in df.index])
    keep_yrs = []
    # find all years that appear >= 365 times
    for yr in yrs:
        if sum([str(idx).startswith(yr) for idx in df.index]) >= 365:
            keep_yrs.append(yr)
    # subset for those
    out = df.loc[[str(idx)[:4] in keep_yrs for idx in df.index]]
    return out


# fit the harmonic regression to each, then
# get the detrended, predicted values for each site
def fit_harmonic_regression(df, response):
    # gather predictors into an array, including:

    # day of total time series (to detrend)
    t = [1 - (max(df.index)-idx)/(
        max(df.index)-min(df.index)) for idx in df.index]

    # sines and cosines of annual and semiannual components
    sin_ann = [np.sin(n) for n in df['ann']]
    cos_ann = [np.cos(n) for n in df['ann']]
    sin_sem = [np.sin(n) for n in df['sem']]
    cos_sem = [np.cos(n) for n in df['sem']]

    X = np.stack((t, sin_ann, cos_ann, sin_sem, cos_sem)).T

    # and grab the y
    y = np.atleast_2d(np.array(df[response])).T

    # build and fit the regression
    reg = linear_model.LinearRegression().fit(X, y)

    return reg


# predict the detrended values
def predict_detrended_vals(df, mod):

    # sines and cosines of annual and semiannual components
    sin_ann = [np.sin(n) for n in df['ann']]
    cos_ann = [np.cos(n) for n in df['ann']]
    sin_sem = [np.sin(n) for n in df['sem']]
    cos_sem = [np.cos(n) for n in df['sem']]
    X = np.stack((sin_ann, cos_ann, sin_sem, cos_sem)).T

    predicted = mod.intercept_ + np.sum(X * mod.coef_[0][1:], axis=1)
    df['pred'] = predicted

# main fn
def process_site_data(csv_filename, response,
                      filter_start_date=None, filter_end_date=None):
    df = pd.read_csv(csv_filename)

    # make datestamp columns in datetime objects
    df['date'] = pd.to_datetime(df['TIMESTAMP'], format='%Y%m%d')

    # subset GPP vars
    gpp = df.loc[:,[*df.columns[[col.startswith('GPP') for col in
                                         df.columns]]]+['date']]
    gpp.set_index('date', inplace=True)

    # drop missing data
    gpp = gpp[gpp[response] != -9999]

    # filter out dates, if needed
    if filter_start_date is not None and filter_end_date is not None:
        # for now, just drop years after 2006 in the forest dataset,
        # because of the big gap caused by the missing 2007-ish data
        # TODO: COME UP WITH BETTER SOLUTION
        gpp = gpp.loc[filter_start_date:filter_end_date]

    gpp = keep_only_full_years(gpp)

    # add numeric day of year columns,
    # and add annual and semi-annual circular time compononet columns (radians)
    gpp['doy'] = [idx.timetuple().tm_yday for idx in gpp.index]
    gpp['ann'] = 2*np.pi*np.array([*gpp['doy']])/365
    gpp['sem'] = 2*np.pi*np.array([*gpp['doy']])/(365/2)

    # fit the regression, then add the predicted, detrended vals col
    reg = fit_harmonic_regression(gpp, response)

    predict_detrended_vals(gpp, reg)

    return gpp



def plot_dfs(response, **kwargs):
    """
    """
    nplots = (len(kwargs))
    fig = plt.figure()
    axs = [fig.add_subplot(nplots, 1, i) for i in range(1, nplots+1)]
    for item, ax in zip(kwargs.items(), axs):
        name, df = item
        df.loc[:,response].plot(ax=ax)
        df.loc[:,'pred'].plot(ax=ax)
        ax.set_title(name)
        ax.set_ylabel('GPP (DT_VUT_50)')

    # add vertical lines for years
    first_days = [df[df.doy == 1].index for name, df in kwargs.items()]
    for n, ax in enumerate(axs):
        first_days_df = first_days[n]
        for first_day in first_days_df:
            ax.plot([first_day]*2, ax.get_ylim(), ':', color='gray', alpha=0.5,
                    linewidth=0.5)
    plt.show()
    return fig



#####################
#### plot Brazil data
#####################

# primary forest
sant_primary = process_site_data(
    './BR/FLX_BR-Sa1_FLUXNET2015_SUBSET_DD_2002-2011_1-4.csv', response_fluxnet,
    '2001-01-01', '2007-01-01')
# logged forest
sant_logged = process_site_data(
    './BR/FLX_BR-Sa3_FLUXNET2015_SUBSET_DD_2000-2004_1-4.csv', response_fluxnet)

# plot
BR_fig = plot_dfs(response_fluxnet,
                  **{'BR_Santarem_forest': sant_primary,
                     'BR_Santarem_logged': sant_logged})


######################
# plot California data
######################

# wetland
CA_mayberry = process_site_data(
    './US/FLX_US-Myb_FLUXNET2015_SUBSET_DD_2010-2014_2-4.csv', response_fluxnet)
# ranch
CA_vaira = process_site_data(
    './US/FLX_US-Var_FLUXNET2015_SUBSET_DD_2000-2014_1-4.csv', response_fluxnet)
# montane forest
CA_blodgett = process_site_data(
    './US/FLX_US-Blo_FLUXNET2015_SUBSET_DD_1997-2007_1-4.csv', response_fluxnet)

# plot
CA_fig = plot_dfs(response_fluxnet,
                  **{'CA_Mayberry': CA_mayberry.iloc[:365*4,:],
                     'CA_Vaira': CA_vaira.iloc[:365*4,:],
                     'CA_Blodgett': CA_blodgett.iloc[:365*4,:]})



#############################
# plot Malaysia and Indo data
#############################

ID = process_site_data('MY_ID/FLX_ID-Pag_FLUXNET-CH4_DD_2016-2017_1-1.csv',
                       response_ch4)
MY = process_site_data('MY_ID/FLX_MY-MLM_FLUXNET-CH4_DD_2014-2015_1-1.csv',
                       response_ch4)
ID_MY_fig = plot_dfs(response_ch4, **{'ID': ID, 'MY': MY})



#####################
# plot Australia data
#####################
AU_fig = plot_dfs(response_fluxnet,
                  **{f.split('-')[1].split('_')[0]:process_site_data(
                f, response_fluxnet) for f in os.listdir('./AU') if 'DD' in f})





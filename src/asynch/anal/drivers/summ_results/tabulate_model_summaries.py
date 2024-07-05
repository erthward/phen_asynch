import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
from affine import Affine
from collections import Counter as C
import os, re, sys

# local imports
sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc/'))
import phen_helper_fns as phf

# set working directory
data_dir = phf.EXTERNAL_RF_DATA_DIR

for coords_as_covars in ['y', 'n']:
    out_dfs = {}
    out_model_dfs = []
    for import_metric in ['permut', 'SHAP']:
        dfs = []
        model_dfs = []
        # loop over vars and neigh_rads (in km)
        for var in ['NIRv', 'SIF']:
            for neigh_rad in [50, 100, 150]:

                # get CSV of Shapley importance values
                df = pd.read_csv(os.path.join(data_dir,
                    'rf_%s_importance_%sCOORDS_%s_%ikm.csv' % (import_metric,
                                                               coords_as_covars,
                                                               var,
                                                               neigh_rad)))
                # rename columns in permut to match those in SHAP importance
                if import_metric == 'permut':
                    df.columns = ['Variable', 'Importance']
                df.loc[:, 'dataset'] = var
                df.loc[:, 'neigh_radius'] = neigh_rad
                dfs.append(df)

                # get the model R2 and RMSE values as well
                output_filename = os.path.join(data_dir,
                        f'ch3_rf_{var}_{neigh_rad}_{coords_as_covars}.Rout')
                with open(output_filename, 'r') as f:
                    text = f.read()
                mse = float(re.search('(?<=OOB prediction error \(MSE\):).*(?=\n)',
                                text).group())
                r2 = float(re.search('(?<=R squared \(OOB\):).*(?=\n)',
                                     text).group())
                model_df = pd.DataFrame({'metric':['R2',
                                                   'MSE'], 'value':[r2, mse]})
                model_df.loc[:, 'dataset'] = var
                model_df.loc[:, 'neigh_radius'] = neigh_rad
                model_dfs.append(model_df)

        all_df = pd.concat(dfs)
        all_model_df = pd.concat(model_dfs)

        out_df = all_df.set_index(['dataset', 'neigh_radius']).pivot(
                    columns=['Variable']).T.reset_index().set_index('Variable').drop(
                    labels=['level_0'], axis=1)
        # drop index name and MultiIndex column names
        out_df.index.name = ''
        out_df.columns.names = ['', '']

        # drop x and y rows, if model did not use coords as covars and this is
        # the SHAP importance metric (NOTE: already automatically dropped for
        # permutation importance metric)
        if coords_as_covars == 'n' and import_metric == 'SHAP':
            out_df = out_df.drop(labels=['x', 'y'], axis=0)

        out_dfs[import_metric] = out_df

        # finalize the model-summary df
        out_model_df = all_model_df.set_index(['dataset', 'neigh_radius']).pivot(
                    columns=['metric']).T.reset_index().set_index('metric').drop(
                    labels=['level_0'], axis=1)
        # convert MSE to RMSE
        out_model_df.loc['MSE', :] = np.sqrt(out_model_df.loc['MSE', :])
        out_model_df.index = ['RMSE', 'R2']
        out_model_df.index.name = ''
        out_model_df.columns.names = ['', '']


    # write to file
    with pd.ExcelWriter(os.path.join(data_dir,
                    'model_summary_table_%sCOORDS.xlsx' % coords_as_covars),
                        engine='xlsxwriter') as w:
        for import_metric, out_df in out_dfs.items():
            out_df.to_excel(w, sheet_name='%s_importance' % import_metric)
            # grab the workseet and apply conditional colors to each column
            worksheet = w.sheets['%s_importance' % import_metric]
            worksheet.conditional_format('B4:G%i' % (4+len(out_df)-1),
                                         {'type': '2_color_scale',
                                          'min_color': '#faeccf',
                                          'max_color': '#f77307',
                                         })
        out_model_df.to_excel(w, sheet_name='model_summary')



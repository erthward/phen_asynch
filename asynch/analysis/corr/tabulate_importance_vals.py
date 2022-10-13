import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
from affine import Affine
from collections import Counter as C
import os

# set working directory
data_dir = '/media/deth/SLAB/diss/3-phn/corr_data/'

out_dfs = {}
for import_metric in ['permut', 'SHAP']:
    dfs = []
    # loop over vars and neigh_rads (in km)
    for var in ['NIRv', 'SIF']:
        for neigh_rad in [50, 100, 150]:

            # get CSV of Shapley importance values
            df = pd.read_csv(os.path.join(data_dir,
                    'rf_%s_importance_%s_%ikm.csv' % (import_metric,
                                                      var,
                                                      neigh_rad)))
            # rename columns in permut to match those in SHAP importance
            if import_metric == 'permut':
                df.columns = ['Variable', 'Importance']
            df.loc[:, 'dataset'] = var
            df.loc[:, 'neigh_radius'] = neigh_rad
            dfs.append(df)

    all_df = pd.concat(dfs)

    out_df = all_df.set_index(['dataset', 'neigh_radius']).pivot(
                columns=['Variable']).T.reset_index().set_index('Variable').drop(
                labels=['level_0'], axis=1)
    # drop index name and MultiIndex column names
    out_df.index.name = ''
    out_df.columns.names = ['', '']
    out_dfs[import_metric] = out_df

# write to file
with pd.ExcelWriter(os.path.join(data_dir,
                                 'var_importance_table.xlsx'),
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

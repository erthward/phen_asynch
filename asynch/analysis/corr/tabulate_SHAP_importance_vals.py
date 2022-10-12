import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
from affine import Affine
from collections import Counter as C
import os

# set working directory
data_dir = '/media/deth/SLAB/diss/3-phn/corr_data/'

dfs = []
# loop over vars and neigh_rads (in km)
for var in ['NIRv', 'SIF']:
    for neigh_rad in [50, 100, 150]:

        # get CSV of Shapley importance values
        df = pd.read_csv(os.path.join(data_dir,
                                       'rf_SHAP_importance_%s_%ikm.csv' % (var,
                                                                           neigh_rad)))
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

# write to file
with pd.ExcelWriter(os.path.join(data_dir, 'SHAP_importance_table.xlsx'),
                    engine='xlsxwriter') as w:
    out_df.to_excel(w, sheet_name='SHAP_importance')

    # grab the workseet and apply conditional colors to each column
    worksheet = w.sheets['SHAP_importance']
    worksheet.conditional_format('B4:G%i' % (4+len(out_df)-1),
                                 {'type': '2_color_scale',
                                  'min_color': '#faeccf',
                                  'max_color': '#f77307',
                                 })

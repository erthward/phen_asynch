import pandas as pd
import numpy as np

# read in evaluation results
nirv_val = pd.read_csv('./FLUXNET_evaluation_results_NIRv.csv')
sif_val = pd.read_csv('./FLUXNET_evaluation_results_SIF.csv')

# combine into single table
nirv_val = nirv_val[pd.notnull(nirv_val['r2'])]
nirv_val.rename(columns={'r2': 'R2_NIRv'}, inplace=True)
sif_val = sif_val[pd.notnull(sif_val['r2'])]
sif_val.rename(columns={'r2': 'R2_SIF'}, inplace=True)
keep_cols = ['id', 'name', 'lon', 'lat', 'igbp', 'mat', 'map']
val = pd.merge(nirv_val[keep_cols+['R2_NIRv']],
               sif_val[keep_cols+['R2_SIF']],
               how='outer',
               on=keep_cols,
              )
val.rename(columns={'id': 'FLUXNET_ID',
                    'name': 'site_name',
                   },
           inplace=True)

# round R2 values to 3 decimals
val['R2_SIF'] = ['%0.3f' % r2 for r2 in val['R2_SIF']]
val['R2_NIRv'] = ['%0.3f' % r2 for r2 in val['R2_NIRv']]

# print number of evaluation datasets for each RS dataset
print(f'\n\n{np.sum(pd.notnull(val["R2_NIRv"]))} val datasets for NIRv')
print(f'\n\n{np.sum(pd.notnull(val["R2_SIF"]))} val datasets for SIF')
print(f'\n\n{len(val)} total datasets')

# write out
val.to_csv('./TABLE_S1_FLUXNET_evaluation_all_results.csv',
           index=False,
          )

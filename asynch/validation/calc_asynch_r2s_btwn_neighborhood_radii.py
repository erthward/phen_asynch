import rioxarray as rxr
import numpy as np
import pandas as pd
import scipy.stats
import itertools
import os, sys, re

data_dir = '/media/deth/SLAB/diss/3-phn/GEE_outputs/final'

vars = ['NIRv_STRICT', 'SIF_STRICT', 'tmmn', 'tmmx', 'pr', 'def', 'cloud']
neighs = [50, 100, 150]
neigh_comps = [*itertools.combinations(neighs, 2)]

var_col = []
neigh_comp_col = []
r2_col = []

for var in vars:
    for n1, n2 in neigh_comps:
        print(f'\ncalculating R2s for {var}: {n1} km and {n2} km...\n')
        r1 = rxr.open_rasterio(os.path.join(data_dir,
                                            f'{var}_asynch_{n1}km.tif'))[0]
        r2 = rxr.open_rasterio(os.path.join(data_dir,
                                            f'{var}_asynch_{n2}km.tif'))[0]
        vals1 = r1.values.ravel()
        vals2 = r2.values.ravel()
        not_nulls = (pd.notnull(vals1))*(pd.notnull(vals2))
        vals1 = vals1[not_nulls]
        vals2 = vals2[not_nulls]
        r2 = scipy.stats.linregress(vals1, vals2)[2]**2
        assert 0 <= r2 <= 1

        var_col.append(var)
        neigh_comp_col.append(f'{n1}|{n2}')
        r2_col.append(r2)

        del r1, r2

df = pd.DataFrame({'var': var_col,
                   'neigh_comp': neigh_comp_col,
                   'r2': r2_col,
                  })
df.to_csv('asynch_neigh_rad_comp_r2s.csv', index=False)


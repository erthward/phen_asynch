import geopandas as gpd
import pandas as pd
import os
import re

# load global adm1 boundaries
raw_adm1 = gpd.read_file('../../data/bounds/ne_10m_admin_1_states_provinces.shp')

# keep only a few columns
adm1 = raw_adm1.loc[:, ['name', 'adm1_code', 'adm0_a3',
                        'area_sqkm', 'latitude', 'longitude', 'geometry']]

# drop Antarctica
adm1 = adm1[adm1['name'] != 'Antarctica']

# derive level 0 adm bounds
adm0 = adm1.dissolve(by='adm0_a3').reset_index()

# save both to file
adm0.to_file('../../data/bounds/global_adm0.shp', index=False)
adm1.to_file('../../data/bounds/global_adm1.shp', index=False)

# load adm1 boundaries for nations with >2.5e6 km^2 land area
select_adm1 = []
for f in os.listdir('../../data/bounds/'):
    if re.search('^gadm.*json$', f):
        select_adm1.append(gpd.read_file(os.path.join('../../data/bounds/',
                                                      f)))
select_adm1 = pd.concat(select_adm1)

# save to file
select_adm1.to_file('../../data/bounds/select_adm1.shp', index=False)

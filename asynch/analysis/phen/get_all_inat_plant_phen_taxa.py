import pyinaturalist as pynat
import numpy as np
import pandas as pd
import geopandas as gpd
import random
import string
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from shapely.geometry import MultiPoint
import contextily as ctx


# set fixed API params
# (NOTE: for details, see return value of `pynat.get_controlled_terms()`)
term_id = 12                # "Plant phenology"
term_value_id = 13          # "Flowering"
quality_grade = 'research'
native = True
captive = False
per_page = 500

# get actual observations
taxa_dict = {'tid': [],
             'name': [],
             'count': [],
            }
curr_page = 0
taxa_pg_ct = 999999
while curr_page < taxa_pg_ct:
    curr_page += 1
    print(f"\tgetting page {curr_page} of valid taxa...")
    taxa = pynat.get_observation_species_counts(term_id=term_id,
                                                term_value_id=term_value_id,
                                                quality_grade=quality_grade,
                                                native=native,
                                                captive=captive,
                                                per_page=per_page,
                                                page=curr_page,
                                               )
    if curr_page == 1:
        taxa_ct = taxa['total_results']
        taxa_pg_ct = int(np.ceil(taxa_ct/per_page))

    # store results
    for taxon in taxa['results']:
        taxa_dict['tid'].append(taxon['taxon']['id'])
        taxa_dict['name'].append(taxon['taxon']['name'])
        taxa_dict['count'].append(taxon['count'])

taxa_df = pd.DataFrame(taxa_dict)
assert len(np.unique(taxa_df['tid'])) == len(taxa_df)
taxa_df.to_csv('all_inat_plant_phen_taxa.csv',
               index=False,
              )

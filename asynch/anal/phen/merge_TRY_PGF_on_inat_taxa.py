import pandas as pd
from collections import Counter as C
from fuzzywuzzy import fuzz, process

# get iNat taxa
taxa = pd.read_csv('./all_inat_plant_phen_taxa.csv')
taxa.loc[:, 'name_low'] = [n.strip().lower() for n in taxa['name'].values]

# get reclassed TRY PGF data
pgf = pd.read_csv('./TRY_PGF_values_reclassed_filtered.csv')

# drop duplicate rows and summarize and decide PGF
pgf = pgf.drop_duplicates(['AccSpeciesID',
                           'SpeciesName',
                           'AccSpeciesName',
                           'forb',
                           'wood',
                           'vine',
                           'epip',
                           'xero',
                          ],
                          ignore_index=True,
                         )
spp_cts = C(pgf.AccSpeciesID.values)

# drop rows missing species names
pgf = pgf[pd.notnull(pgf['SpeciesName']) & pd.notnull(pgf['AccSpeciesName'])]

# NOTE: only 3.26% of rows in this original PGF table have a True for both
#       'forb' and 'wood', and those 2 PGFs cover the vast majority of the taxa
#       we'll by analyzing, so rather than worry about reconciling PGF now
#       we'll just leave it tacked onto the taxa table and then worry about
#       reconciling that if/as necessary at the end of the analysis

# reconcile rows with 2 different species names
# (defaulting the AccSpeciesName value)
name_to_use = []
for i, row in pgf.iterrows():
    accname = row['AccSpeciesName'].strip().lower()
    name = row['SpeciesName'].strip().lower()
    names_equal = accname == name
    if names_equal or accname in taxa['name_low'].values:
        name_to_use.append(accname)
    elif name in taxa['name_low'].values:
        name_to_use.append(name)
    else:
        name_to_use.append(accname)
pgf.loc[:, 'name_low'] = name_to_use

# merge
merged = pd.merge(taxa, pgf, on='name_low', how='left')

# summarize PGFs by taxon
pgf_cols = ['forb', 'wood', 'vine', 'epip', 'xero']
pgf_strs = []
for tid in taxa['tid'].values:
    pgfs = merged[merged['tid'] == tid].loc[:, pgf_cols]
    pgf_str = ''.join([s[0] for s,
                              v in pgfs.sum(axis=0).to_dict().items() if v !=0])
    pgf_strs.append(pgf_str)
taxa.loc[:, 'pgf'] = pgf_strs

# write to disk
taxa.to_csv('all_inat_plant_phen_taxa_w_TRY_pgf.csv',
            index=False,
           )

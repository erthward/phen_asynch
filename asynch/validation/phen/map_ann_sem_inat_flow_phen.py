import pyinaturalist as pynat
import numpy as np
import pandas as pd
import geopandas as gpd
import random
import string
import time
import os
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from shapely.geometry import MultiPoint
import contextily as ctx


# TODO:

    # add numerous QC filters
    # (read through and understand all possible args, using that construct
    #  robust filtering)

    # how to choose taxa? or just take everything?

# read in all taxa with plant phenology data
taxa = pd.read_csv('./all_inat_plant_phen_taxa.csv')
assert len(np.unique(taxa['tid'])) == len(taxa)
print(f"\n{len(taxa)} total taxa with flowering phenology info in iNat\n\n")

# compare and subset against already-processed taxa
processed_file_writemode = 'w'
processed_taxa_filename = 'inat_flower_phen_res.json'
if os.path.isfile(processed_taxa_filename):
    processed_file_writemode = 'a'
    print('\treading already-processed taxa...\n')
    processed_taxa = gpd.read_file(processed_taxa_filename)
    taxa = taxa[~taxa['tid'].isin(processed_taxa['tid'])]
print(f"\n{len(taxa)} taxa remaining to be processed\n\n")

# what is the max number of points we want to pull per species?
max_points_per_species = 2000

# set fixed API params
# (NOTE: for details, see return value of `pynat.get_controlled_terms()`)
term_id = 12                # "Plant phenology"
term_value_id = 13          # 13: "Flowering", 14: "Fruting",
                            # 15: "Flower Budding",
                            # 21: "No Evidence of Flowering"
date_file = 'observed'
interval = 'week_of_year'
quality_grade = 'research'
per_page = 200
positional_accuracy = 1000
native = True

# regression df
reg_df = pd.DataFrame.from_dict({'woy': [*range(1,54)]})
reg_df['woy_sin_ann'] = np.sin(reg_df['woy']/53*2*np.pi)
reg_df['woy_cos_ann'] = np.cos(reg_df['woy']/53*2*np.pi)
reg_df['woy_sin_sem'] = np.sin(reg_df['woy']/53*2*2*np.pi)
reg_df['woy_cos_sem'] = np.cos(reg_df['woy']/53*2*2*np.pi)

# results data struct
res_dict ={'tid': [],
           'name': [],
           'count': [],
           'r2_ann': [],
           'r2_sem': [],
           'geom': [],
          }

# simple function to scale array between 0 and n
def scale_arr(arr, max_val):
    return (arr - np.min(arr))/(np.max(arr) - np.min(arr)) * max_val


# simple function to get a random hex color (not too light!)
def get_hex_color():
    hex_choices = random.choices(string.hexdigits[3:], k=6)
    return f"#{''.join(hex_choices)}".upper()


# for each taxon, get histogram
failed = False
colors = []
for i, row in taxa.iterrows():
    success = False
    trys = 0
    while not success and not failed:
        try:
            time.sleep(1.5)
            tid = row['tid']
            tax_name = row['name']
            max_count = row['count']
            print(f"\n\n{'.'*80}\nprocessing {tax_name}...")

            # calculate area and centroid of geographic range of observations

            # get observation histogram
            hist = pynat.get_observation_histogram(taxon_id=tid,
                                                   term_id=term_id,
                                                   term_value_id=term_value_id,
                                                   date_file=date_file,
                                                   interval=interval,
                                                   quality_grade=quality_grade,
                                                   native=native,
                                                  )

            # compare R2s between annual and semi-annual harmonic regressions
            reg_df[tid] = [*hist.values()]
            reg_ann = LinearRegression().fit(X=reg_df.loc[:, ['woy_sin_ann', 'woy_cos_ann']],
                                             y=reg_df[tid].values,
                                            )
            r2_ann = reg_ann.score(X=reg_df.loc[:, ['woy_sin_ann', 'woy_cos_ann']],
                                   y=reg_df[tid].values,
                                  )
            reg_sem = LinearRegression().fit(X=reg_df.loc[:, ['woy_sin_sem', 'woy_cos_sem']],
                                             y=reg_df[tid].values,
                                            )
            r2_sem = reg_sem.score(X=reg_df.loc[:, ['woy_sin_sem', 'woy_cos_sem']],
                                   y=reg_df[tid].values,
                                  )

            # get actual observations
            curr_page = 0
            locs = []
            obs_pg_ct = 999999
            while curr_page < obs_pg_ct:
                curr_page += 1
                print(f"\tgetting page {curr_page} of observation locations...")
                obs = pynat.get_observations(taxon_id=tid,
                                             term_id=term_id,
                                             term_value_id=term_value_id,
                                             date_file=date_file,
                                             quality_grade=quality_grade,
                                             per_page=per_page,
                                             page=curr_page,
                                             native=native,
                                            )

                # NOTE: coords returned as lat,lon instead of lon,lat!
                locs.extend([o['location'][::-1] for o in obs['results'] if
                             (o['positional_accuracy'] is not None and
                              o['positional_accuracy']<=positional_accuracy)])
                if curr_page == 1:
                    # NOTE: taking a maximum number of points per species
                    #       to speed things up, at least for now
                    obs_ct = np.min((obs['total_results'], max_points_per_species))
                    obs_pg_ct = int(np.ceil(obs_ct/per_page))

            # store results
            res_dict['tid'].append(tid)
            res_dict['name'].append(tax_name)
            res_dict['count'].append(len(locs))
            res_dict['r2_ann'].append(r2_ann)
            res_dict['r2_sem'].append(r2_sem)
            res_dict['geom'].append(MultiPoint(locs))

            # plot hist and both regressions
            max_r2 = np.max([r2_ann, r2_sem])
            line_types = {freq: [':', '-'][int(r2 == max_r2)] for
                                    freq, r2 in zip(['ann', 'sem'], [r2_ann, r2_sem])}
            color = get_hex_color()
            fig, ax = plt.subplots(1)
            ax.plot([*range(len(hist))], [*hist.values()], ':k', label='observations')
            pred_ann = np.sum(reg_df.loc[:, ['woy_sin_ann',
                                             'woy_cos_ann']] * reg_ann.coef_, axis=1)
            pred_sem = np.sum(reg_df.loc[:, ['woy_sin_sem',
                                             'woy_cos_sem']] * reg_ann.coef_, axis=1)
            ax.plot([*range(len(hist))],
                    scale_arr(pred_ann, np.max([*hist.values()])),
                    line_types['ann'],
                    color=color,
                    label='annual fit',
                   )
            ax.plot([*range(len(hist))],
                    scale_arr(pred_sem, np.max([*hist.values()])),
                    line_types['sem'],
                    color=color,
                    label='semiannual fit',
                   )
            fig.legend()
            ax.set_xlabel('week of year')
            ax.set_ylabel('number of observations')
            ax.set_title((f"TID {tid}: {tax_name}\n"
                f"({np.sum([*hist.values()])} total flowering observations)"))
            fig.savefig(f"./plots/TID_{tid}_{tax_name.replace(' ', '_')}.png",
                        dpi=300,
                       )
            colors.append(color)
            success = True
        except Exception as e:
            trys += 1
            if trys > 3:
                print('\ntoo many trys; breaking loop...')
                failed = True
                break
            print(f"EXCEPTION THROWN: {e}")
            print('\n\nretrying...')


res_df = pd.DataFrame.from_dict(res_dict)
res_gdf = gpd.GeoDataFrame(res_df,
                           crs=4326,
                           geometry='geom',
                          )
assert len(np.unique(res_gdf['tid'])) == len(res_gdf)
res_gdf['r2_ratio'] = res_gdf['r2_sem']/res_gdf['r2_ann']
if failed:
    res_gdf.iloc[:-1, :].to_file(processed_taxa_filename,
                                 mode=processed_file_writemode)
else:
    res_gdf.to_file(processed_taxa_filename,
                    mode=processed_file_writemode)
    fig, ax = plt.subplots(1)
    res_gdf['color'] = colors
    res_gdf.to_crs(3857).plot(alpha=0.7,
                              color=res_gdf['color'],
                              edgecolor='black',
                              markersize=10,
                              ax=ax,
                              legend='name',
                             )
    ctx.add_basemap(ax=ax)
    fig.savefig("inat_flower_phen_obs_locs.png",
            dpi=600,
           )


# do places with higher asynchrony have a higher percentage of semi-annual
# flowering phenologies?

# how does that percentage map/krig out?

# TODO: mapping/analysis thoughts:
    # - calculate ratio of r2_sem/r2_ann
    # - plot distribution of all centroids, then overplot all centroids with ratio>1
    # - regression of ratio and asynch values to abs(lat)? point asynch? etc


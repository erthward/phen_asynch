import pyinaturalist as pynat
import numpy as np
import pandas as pd
import geopandas as gpd
import alphashape
from shapely.geometry import MultiPolygon
import matplotlib.pyplot as plt
import datetime
import time

# which taxon ID to use?
# Dipterostemon capitatus: 2nd most samples of any taxon (>=11963),
# and has a fairly complex range
#tid_to_use = 1196784
# Geranium robertianum: widespread in Europe and N. Africa
tid_to_use = 55925

# read in taxa table
taxa = pd.read_csv('./all_inat_plant_phen_taxa_w_TRY_pgf.csv')
taxon = taxa[taxa['tid'] == tid_to_use]['name'].values[0]

# alpha-hull alpha param value to use?
alpha = 0.25 # this is the mid-range value we used in Fig. 4 analysis

# get all observations
curr_page = 0
coords = []
obs_dates = []
per_page = 200
obs_pg_ct = 999999
term_id = 12
term_value_id = 13
date_file = 'observed'
quality_grade = 'research'
native = True
min_pos_acc = 1000
while curr_page < obs_pg_ct:
    curr_page += 1
    success = False
    while not success:
        try:
            print(f"\tgetting page {curr_page} of observation locations...")
            obs = pynat.get_observations(taxon_id=tid_to_use,
                                         term_id=term_id,
                                         term_value_id=term_value_id,
                                         date_file=date_file,
                                         quality_grade=quality_grade,
                                         per_page=per_page,
                                         page=curr_page,
                                         native=native,
                                        )
            success = True
        except Exception as e:
            print(f"ERROR THROWN: {e}\n\n\nwaiting 1 minute...\n")
            time.sleep(60)
    # only keep points with adequate positional accuracy
    for o in obs['results']:
        if (o['positional_accuracy'] is not None and
            o['positional_accuracy'] <= min_pos_acc):
            # get date of observation
            obs_date = o['observed_on']
            if isinstance(obs_date, str):
                obs_date = datetime.datetime.strptime(obs_date, '%Y-%m-%d')
            obs_dates.append(obs_date)
            # get coords
            # NOTE: returned as lat,lon instead of lon,lat!
            coord = o['location'][::-1]
            coords.append(coord)
    if curr_page == 1:
        # NOTE: taking a maximum number of points per species
        #       to speed things up, at least for now
        obs_ct = obs['total_results']
        # NOTE: for some reason we only seem to be able to get up through pg 50
        obs_pg_ct = np.min((int(np.ceil(obs_ct/per_page)), 50))

assert len(coords) == len(obs_dates)

# loop over n values, subsetting first n observations and fitting alpha hulls
# NOTE: only got about 7363 coordinates because of positional accuracy
hulls = []
ns = [500, 1000, 2000, 3000, 4000, len(coords)]
coords = np.array(coords)
for n in ns:
    alpha_hull = alphashape.alphashape(coords[:n,:],
                                       alpha=alpha,
                                      )
    if isinstance(alpha_hull, MultiPolygon):
        pass
    else:
        alpha_hull = MultiPolygon([alpha_hull])
    hulls.append(alpha_hull)

# map them all
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
fig = plt.figure(figsize=(10,10))
axs = fig.subplots(2,3)
for n, ax, hull in zip(ns, axs.ravel(), hulls):
    ax.set_title(f"first {n} points", fontdict={'fontsize': 14})
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(())
    ax.set_yticks(())
    world.plot(color='none', linewidth=0.5, zorder=0, ax=ax)
    ax.scatter(coords[:n, 0], coords[:n, 1], marker='.', s=1, zorder=1)
    for poly in hull:
        x, y = poly.exterior.xy
        ax.plot(x, y, color='orange', alpha=0.9, linewidth=1, zorder=2)
    ax.set_xlim(1.1 * np.min(coords[:,0]), 0.9 * np.max(coords[:, 0]))
    ax.set_ylim(0.9 * np.min(coords[:,1]), 1.1 * np.max(coords[:, 1]))
fig.suptitle(f"{taxon}: TID {tid_to_use}", size=19)
fig.savefig((f"TID_{tid_to_use}_{taxon.replace(' ', '_')}"
             f"_hulls_vs_n_pts_ALPHA{str(alpha)}.png"),
            dpi=600,
           )


# write observation coords to CSV
df = pd.DataFrame({'lat': coords[:, 0],
                   'lon': coords[:, 1],
                   'date': obs_dates,
                  })
df.to_csv(f"TID_{tid_to_use}_{taxon.replace(' ', '_')}_all_flow_phen_pts.csv",
          index=False,
         )



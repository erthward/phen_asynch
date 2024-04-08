from shapely.geometry import Point, MultiPolygon
from sklearn.linear_model import LinearRegression
from pyinaturalist import clear_cache
import matplotlib.pyplot as plt
import pyinaturalist as pynat
import contextily as ctx
import geopandas as gpd
import pandas as pd
import numpy as np
import alphashape
import subprocess
import datetime
import diptest
import random
import string
import time
import os

# silence the unhelpful shapely deprecation warnings
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning) 


#------------------------------------------------------------------------------
# path mgmt
#------------------------------------------------------------------------------
data_dir = '/media/deth/SLAB/diss/3-phn/inat/'


#------------------------------------------------------------------------------
# data load
#------------------------------------------------------------------------------

# set numpy seed
np.random.seed(1)

# read in all taxa with plant phenology data
taxa = pd.read_csv('./all_inat_plant_phen_taxa_w_TRY_pgf.csv')
assert len(np.unique(taxa['tid'])) == len(taxa)
print('\n\n')
print(f"\n{len(taxa)} total taxa with flowering phenology info in iNat\n")

# subset to only taxa with the minimum requisite number of observations
min_obs = 50
taxa = taxa[taxa['count'] >= min_obs]
print(f"\n{len(taxa)} taxa with >= {min_obs} flowering phenology observations\n")

# compare and subset against already-processed taxa
processed_file_writemode = 'w'
processed_taxa_filename = 'inat_flower_phen_results.json'
if os.path.isfile(processed_taxa_filename):
    processed_file_writemode = 'a'
    print('\treading already-processed taxa...')
    processed_taxa = gpd.read_file(processed_taxa_filename)
    taxa = taxa[~taxa['tid'].isin(processed_taxa['tid'])]
# reset the index (so that iterrows gives incrementing ints starting from 0)
taxa = taxa.reset_index()
print((f"\n{len(taxa)} taxa remaining to be processed "
       f"({len(processed_taxa)} already complete)\n\n"))

# regression df
# NOTE: has to be 53 weeks, to match iNat output
#       (since a week technically has 52.143 weeks on non-leap years)
reg_df = pd.DataFrame.from_dict({'woy': [*range(1,54)]})
reg_df['woy_sin_ann'] = np.sin(reg_df['woy']/53*2*np.pi)
reg_df['woy_cos_ann'] = np.cos(reg_df['woy']/53*2*np.pi)
reg_df['woy_sin_sem'] = np.sin(reg_df['woy']/53*2*2*np.pi)
reg_df['woy_cos_sem'] = np.cos(reg_df['woy']/53*2*2*np.pi)

# construct canonical 1x and 2x sinusoidal curves, against which data will be
# compared using earth-mover's distance
minmax_scale = lambda a: (a-np.min(a))/(np.max(a)-np.min(a))
sin_1x = minmax_scale(np.cos(np.linspace(np.pi, 3*np.pi, 53)))
sin_2x = minmax_scale(np.cos(np.linspace(np.pi, (3+2)*np.pi, 53)))


#------------------------------------------------------------------------------
# data-collection and processing params
#------------------------------------------------------------------------------
# max positional accuracy value allowed
# NOTE: accuracy values are expressed as "accurate to within <VAL> meters"
#       (from iNat site: "distance in meters that includes the entire area
#        where the observation could have taken place")
max_pos_acc_val = 1000

# max number of points we want to use to construct a taxon's alpha hull
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
native = True
captive = False

# sigma on the normal dist used to add noise to discrete data for diptest input
diptest_sigma = 1

# alpha value to use for the alpha-hull algo
# NOTE: 0.75 is the mid value of the 3 values used in the climate-distance
#       asynchrony analysis (Fig. 4)
alpha = 0.75


#------------------------------------------------------------------------------
# helper fns
#------------------------------------------------------------------------------

def make_obs_dict():
    '''
    make a dict to store observation-level data for a single taxon
    '''
    obs_dict = {'datetime': [],
                'doy':      [],
                'doy_circ': [],
                'geometry': [],
               }
    return obs_dict


def scale_arr(arr, max_val):
    '''
    scale array between 0 and n
    '''
    return (arr - np.min(arr))/(np.max(arr) - np.min(arr)) * max_val


def get_hex_color():
    '''
    get a random hex color (not too light!)
    '''
    hex_choices = random.choices(string.hexdigits[3:], k=6)
    return f"#{''.join(hex_choices)}".upper()


def calc_r2_adj(r2, reg, n):
    '''
    calculate the adjusted R^2 for the given R^2, regression, and n
    (where reg is a sklearn.LinearRegression object)
    '''
    p = reg.coef_.size + int(reg.intercept_ != 0)
    r2_adj = 1 - ((1 - r2) * ((n-1)/(n-p-1)))
    return r2_adj


def calc_euc_dist(a1, a2):
    '''
    calc the n-dimensional Euclidean distance between two n-dimensional arrays
    '''
    return np.sqrt(np.sum([(v1-v2)**2 for v1, v2 in zip(a1, a2)]))


def rotate_time_series_to_min(ts,
                              nsteps=7,
                             ):
    '''
    shift time series so that it starts on its min value
    (or which of its min values is just before the greatest increase over the
    next n steps, if it has more than 1; e.g., zero-inflated observation data)
    '''
    if not isinstance(ts, np.ndarray):
        ts = np.array([*ts])
    minval = np.min(ts)
    #print('ts: ', ts)
    minidx = np.where(ts == minval)[0]
    #print('minidx: ', minidx)
    if len(minidx) == 1:
        #print('just 1')
        cutidx = minidx[0]
    else:
        #print('>1')
        tscat=  np.concatenate([ts]*2)
        minidx_nsteps_slope = (tscat[minidx+nsteps] - tscat[minidx])/nsteps
        #print('minidx_nsteps_slope: ', minidx_nsteps_slope)
        # NOTE: if there is more than one spot with the maximum observed
        #       n-step slope right after a min value then this will just default
        #       to the first one, which should be good enough
        minidx_maxdif = np.argmax(minidx_nsteps_slope)
        #print('minidx_maxdif: ', minidx_maxdif)
        cutidx = minidx[minidx_maxdif]
    #print('cutidx: ', cutidx)
    ts_rot = np.concatenate((ts[cutidx:], ts[:cutidx]))
    return ts_rot


def run_diptest_in_R(vals, is_hist=True, noise_sigma=0):
    '''
    use the histogram values to run the dip test in R

    If input vals are samples from the distribution to be tested then set
    is_hist to False.

    If input vals are densities within the subsequent bins of a histogram
    calculated from samples from that distribution then set is_hist to True.
    '''
    assert not np.all(np.array(vals)==0), 'cannot run dip test on all zeros!'
    tmp_filename = 'diptest_data.tmp'
    if is_hist:
        # rotate
        hist_rot = rotate_time_series_to_min(vals)
        # turn into hist samples
        samp = [i for i,v in enumerate(hist_rot) for _ in range(v)]
        # add noise
        samp += np.random.normal(0, noise_sigma, len(samp))
    else:
        samp = vals
    # cast as data.frame
    df = pd.DataFrame.from_dict({'samp': samp})
    # save to 'diptest_data.tmp'
    df.to_csv(tmp_filename, index=False)
    # call Rscript and get results as dict
    out = subprocess.getoutput('Rscript --vanilla run_diptest.r')
    stats = out.split('\n')
    res = {s.split(': ')[0]: float(s.split(': ')[1]) for s in stats}
    # delete tmp file
    os.remove(tmp_filename)
    # return results
    return res['dip'], res['p']



#------------------------------------------------------------------------------
# get histograms and observations
#------------------------------------------------------------------------------
# overall results data struct 
res_dict ={'tid': [],
           'name': [],
           'hist_count': [],
           'hist_count_expec': [],
           'obs_count': [],
           'r2_ann': [],
           'r2_sem': [],
           'dip_stat': [],
           'dip_pval': [],
           'dist_1xsin': [],
           'dist_2xsin': [],
           'color': [],
           'geometry': [],
          }

# for each taxon, get histogram
total_failure = False
runtimes = []
for i, row in taxa.iterrows():
    obs_dict = make_obs_dict()
    try:
        start_time = time.time()
        tid = row['tid']
        tax_name = row['name']
        print(f"\n\n{'.'*80}\nprocessing {tax_name} ({i+1} of {len(taxa)})...")

        # get observation histogram
        hist_success = False
        n_trys=0
        while not hist_success:
            try:
                hist = pynat.get_observation_histogram(taxon_id=tid,
                                                       term_id=term_id,
                                                       term_value_id=term_value_id,
                                                       date_file=date_file,
                                                       interval=interval,
                                                       quality_grade=quality_grade,
                                                       native=native,
                                                       captive=captive,
                                              )
                hist_success = True
                print((f"\n\t{np.sum([*hist.values()])} observations returned "
                       f"(vs. {taxa[taxa['tid'] == tid]['count'].values[0]} "
                        "expected based on initial taxa table).\n"))
            except Exception as e:
                print((f"Error thrown during histogram API call: {e}\n\n"
                       f"Waiting {60 * n_trys} sec, then trying again...\n"))
                time.sleep(60 * n_trys)
                n_trys += 1

        # if hist is all zeros then skip and return missing data
        hist_vals = [*hist.values()]
        if np.sum(hist_vals) == 0:
            res_dict['tid'].append(tid)
            res_dict['name'].append(tax_name)
            res_dict['hist_count'].append(np.nan)
            res_dict['hist_count_expec'].append(row['count'])
            res_dict['obs_count'].append(np.nan)
            res_dict['r2_ann'].append(np.nan)
            res_dict['r2_sem'].append(np.nan)
            res_dict['dip_stat'].append(np.nan)
            res_dict['dip_pval'].append(np.nan)
            res_dict['dist_1xsin'].append(np.nan)
            res_dict['dist_2xsin'].append(np.nan)
            res_dict['color'].append(np.nan)
            res_dict['geometry'].append(MultiPolygon())

            # store and save taxon observations
            obs_dict['datetime'].extend([pd.NaT])
            obs_dict['doy'].extend([np.nan])
            obs_dict['doy_circ'].extend([np.nan])
            obs_dict['geometry'].extend([Point()])

        # get actual observations
        else:
            obs_success = False
            n_trys = 0
            while not obs_success:
                try:
                    curr_page = 0
                    coords = []
                    dates = []
                    doys = []
                    doys_circ = []
                    obs_pg_ct = 999999
                    while curr_page < obs_pg_ct:
                        curr_page += 1
                        if curr_page == 1:
                            print("\tgetting observations, page 1...")
                        else:
                            print((f"\tgetting observations, page {curr_page} "
                                   f"of {obs_pg_ct}..."))
                        obs = pynat.get_observations(taxon_id=tid,
                                                     term_id=term_id,
                                                     term_value_id=term_value_id,
                                                     date_file=date_file,
                                                     quality_grade=quality_grade,
                                                     per_page=per_page,
                                                     page=curr_page,
                                                     native=native,
                                                     captive=captive,
                                                    )
                        for o in obs['results']:
                            pos_acc_ok = (o['positional_accuracy'] is not None and
                                          o['positional_accuracy']<=max_pos_acc_val)
                            if pos_acc_ok:
                                # get coords
                                # NOTE: given in lat,lon order rather than lon,lat!
                                obs_coords = o['location'][::-1]
                                coords.append(obs_coords)
                                # get datetime
                                obs_date = o['observed_on']
                                if isinstance(obs_date, str):
                                    obs_date = datetime.datetime.strptime(obs_date,
                                                                        '%Y-%m-%d')
                                assert isinstance(obs_date, datetime.datetime)
                                dates.append(obs_date)
                                # convert to doy in rads
                                obs_yr = obs_date.year
                                days_in_yr = datetime.datetime.strptime(
                                                    f"{obs_yr}-12-31",
                                                    '%Y-%m-%d').timetuple().tm_yday
                                doy = obs_date.timetuple().tm_yday
                                doy_circ = (doy/days_in_yr) * 2 * np.pi
                                doys.append(doy)
                                doys_circ.append(doy_circ)
                        if curr_page == 1:
                            # NOTE: taking a maximum number of points per species
                            #       to speed things up, at least for now
                            obs_ct = np.min((obs['total_results'],
                                             max_points_per_species))
                            obs_pg_ct = int(np.ceil(obs_ct/per_page))
                    assert len(coords) == len(dates) == len(doys) == len(doys_circ)
                    obs_success = True
                except Exception as e:
                    print((f"Error thrown during observations API call: {e}\n\n"
                           f"Waiting {60 * n_trys} sec, then trying again...\n"))
                    time.sleep(60 * n_trys)
                    n_trys += 1

            # compare R2s between annual and semi-annual harmonic regressions
            reg_ann = LinearRegression().fit(X=reg_df.loc[:, ['woy_sin_ann', 'woy_cos_ann']],
                                             y=hist_vals,
                                            )
            r2_ann = reg_ann.score(X=reg_df.loc[:, ['woy_sin_ann', 'woy_cos_ann']],
                                   y=hist_vals,
                                  )
            r2_ann_adj = calc_r2_adj(r2_ann, reg_ann, len(reg_df))
            reg_sem = LinearRegression().fit(X=reg_df.loc[:, ['woy_sin_sem', 'woy_cos_sem']],
                                             y=hist_vals,
                                            )
            r2_sem = reg_sem.score(X=reg_df.loc[:, ['woy_sin_sem', 'woy_cos_sem']],
                                   y=hist_vals,
                                  )
            r2_sem_adj = calc_r2_adj(r2_sem, reg_sem, len(reg_df))
            print(f"\n\n\tR^2 ratio: {r2_sem/r2_ann}")

            # get the stat and P-value for a dip test on the histogram
            # NOTE: dip test appears to fail on int (i.e., discrete) data,
            #       but week numbers are coarse approximations of the true
            #       continuous-time timepoint at which each observation was made,
            #       so adding some noise is totally reasonable and if anything should
            #       only make our result more conservative on histograms with some true
            #       indication of multi-modality
            dip, pval = run_diptest_in_R(hist_vals, noise_sigma=diptest_sigma)
            print(f'\tdip: {np.round(dip, 2)} (P-value: {np.round(pval, 2)})\n')

            # get the Euclidean distance between the histogram (rotated to its min
            # value before the greatest slope over the 7 following weeks) and both
            # 1x and 2x sine curves
            eud1 = calc_euc_dist(sin_1x,
                                 minmax_scale(rotate_time_series_to_min(hist_vals)))
            eud2 = calc_euc_dist(sin_2x,
                                 minmax_scale(rotate_time_series_to_min(hist_vals)))

            # calculate alpha of observation coordinates
            hull =  alphashape.alphashape(np.array(coords), alpha=alpha)
            if not isinstance(hull, MultiPolygon):
                hull = MultiPolygon([hull])
            assert isinstance(hull, MultiPolygon)

            # assign a random color to this taxon
            color = get_hex_color()

            # plot hist and both regressions
            max_r2 = np.max([r2_ann, r2_sem])
            line_types = {freq: [':', '-'][int(r2 == max_r2)] for
                                    freq, r2 in zip(['ann', 'sem'], [r2_ann, r2_sem])}
            fig = plt.figure(figsize=(4,4))
            ax = fig.add_subplot(1,1,1)
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
            fig_filename = os.path.join(data_dir,
                        f"phen_plots/TID_{tid}_{tax_name.replace(' ', '_')}.png")
            fig.savefig(fig_filename, dpi=400)
            plt.close('all')

            # store results
            res_dict['tid'].append(tid)
            res_dict['name'].append(tax_name)
            res_dict['hist_count'].append(np.sum(hist_vals))
            res_dict['hist_count_expec'].append(row['count'])
            res_dict['obs_count'].append(len(coords))
            res_dict['r2_ann'].append(r2_ann)
            res_dict['r2_sem'].append(r2_sem)
            res_dict['dip_stat'].append(dip)
            res_dict['dip_pval'].append(pval)
            res_dict['dist_1xsin'].append(eud1)
            res_dict['dist_2xsin'].append(eud2)
            res_dict['color'].append(color)
            res_dict['geometry'].append(hull)

            # store taxon observations
            obs_dict['datetime'].extend([str(d).replace(' ', 'T') for d in dates])
            obs_dict['doy'].extend(doys)
            obs_dict['doy_circ'].extend(doys_circ)
            obs_dict['geometry'].extend([Point(c) for c in coords])

        # save all data
        obs_df = pd.DataFrame.from_dict(obs_dict)
        obs_gdf = gpd.GeoDataFrame(obs_df,
                                   geometry='geometry',
                                   crs=4326,
                                  )
        obs_filename = os.path.join(data_dir,
                        f"obs_data/TID_{tid}_{tax_name.replace(' ', '_')}.json")
        obs_gdf.to_file(obs_filename)

        # save what data we have thus far
        res_df = pd.DataFrame.from_dict(res_dict)
        res_gdf = gpd.GeoDataFrame(res_df,
                                   geometry='geometry',
                                   crs=4326,
                                  )
        assert len(np.unique(res_gdf['tid'])) == len(res_gdf)
        res_gdf['r2_ratio'] = res_gdf['r2_sem']/res_gdf['r2_ann']
        res_gdf_curr_it = res_gdf[res_gdf['name'] == tax_name]
        assert len(res_gdf_curr_it) == 1
        assert res_gdf_curr_it['tid'].values[0] == tid
        res_gdf_curr_it.to_file(processed_taxa_filename,
                                mode=processed_file_writemode)

        # clear cache every 10 calls, so laptop memory doesn't top out
        if i % 10 == 0:
            print("\n\n\tCLEARING CACHE...\n\n")
            clear_cache()

        stop_time = time.time()
        time_diff = stop_time - start_time
        runtimes.append(time_diff)
        print(f"\n\trolling-average runtime: {np.round(np.sum(runtimes)/(i+1), 1)} sec per taxon\n")

    except Exception as e:
        print(f"EXCEPTION THROWN: {e}")
        total_failure = True
        break


#------------------------------------------------------------------------------
# analyze results
#------------------------------------------------------------------------------
if total_failure:
    print("\n\nTOTAL FAILURE! Interim results saved. Debug and rerun.\n\n")
else:
    fig, ax = plt.subplots(1)
    res_gdf.to_crs(3857).plot(alpha=0.7,
                              color=res_gdf['color'],
                              edgecolor='black',
                              markersize=10,
                              ax=ax,
                              legend='name',
                             )
    try:
        ctx.add_basemap(ax=ax)
    except Exception as e:
        print("\n\n\tCONTEXTILY ERROR; NO BASEMAP\n\n")
        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        world = world[world['continent'] != 'Antarctica']
        world.to_crs(3857).plot(color='none', edgecolor='black', linewidth=0.5,
                                alpha=0.5, ax=ax, zorder=0)
    fig.savefig("inat_flower_phen_obs_locs.png", dpi=600)


# do places with higher asynchrony have a higher percentage of semi-annual
# flowering phenologies?

# how does that percentage map/krig out?

# TODO: mapping/analysis thoughts:
    # - calculate ratio of r2_sem/r2_ann
    # - plot distribution of all centroids, then overplot all centroids with ratio>1
    # - regression of ratio and asynch values to abs(lat)? point asynch? etc





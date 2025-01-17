from shapely.geometry import Point, MultiPolygon
from sklearn.neighbors import KernelDensity
from pyinaturalist import clear_cache
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import pyinaturalist as pynat
import contextily as ctx
import geopandas as gpd
import pandas as pd
import numpy as np
import xyzservices
import alphashape
import subprocess
import datetime
import random
import string
import time
import sys
import os

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/src/etc'))
import phen_helper_fns as phf

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
taxa = pd.read_csv('./all_inat_plant_phen_taxa.csv')
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
    print((f"\n{len(taxa)} taxa remaining to be processed "
           f"({len(processed_taxa)} already complete)\n\n"))
    # reset index, so that iterrows gives incrementing ints starting from 0
    taxa = taxa.reset_index()


#------------------------------------------------------------------------------
# data-collection and processing params
#------------------------------------------------------------------------------
# max positional accuracy value allowed
# NOTE: accuracy values are expressed as "accurate to within <VAL> meters"
#       (from iNat site: "distance in meters that includes the entire area
#        where the observation could have taken place")
max_pos_acc_val = 1000

# max number of points we want to use to construct a taxon's alpha hull
max_points_per_species = 5000

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

# alpha value to use for the alpha-hull algo
# NOTE: 0.75 is the mid value of the 3 values used in the climate-distance
#       asynchrony analysis (Fig. 4)
alpha = 0.75

# KDE bandwidth, min peak height, and number of permutations for procedure
# estimating number of histogram peaks and its p-value
npeak_bandwidth = 5
npeak_minheight = 0.6
npeak_nperm = 100

# whether to re-plot the histogram plots
replot_hists = True

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


def get_hex_color():
    '''
    get a random hex color (not too light!)
    '''
    hex_choices = random.choices(string.hexdigits[3:], k=6)
    return f"#{''.join(hex_choices)}".upper()


def rotate_time_series_to_min(ts,
                              nsteps=7,
                              return_cutidx=False,
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
    if return_cutidx:
        return ts_rot, cutidx
    else:
        return ts_rot


def calc_n_kde_peaks(hist,
                     kde_bandwidth=5,
                     peak_minheight=0.6,
                     plot=False,
                     ax=None,
                    ):
    '''
    combine a few metrics to estimate the number of peaks in a histogram
    '''
    # draw samples representing the histogram
    samps = [i for i, v in enumerate(hist) for _ in range(v)]
    S = np.array(samps).reshape((-1, 1))
    # calculate original KDE
    kde = KernelDensity(bandwidth=kde_bandwidth).fit(S)
    # get fitted density
    xs = np.arange(0, len(hist)+0.1, 0.1).reshape((-1, 1))
    log_dens = kde.score_samples(xs)
    dens = np.exp(log_dens)
    minmax_dens = (dens - np.min(dens))/(np.max(dens) - np.min(dens))
    if plot:
        minmax_hist = (hist - np.min(hist))/(np.max(hist) - np.min(hist))
        if ax is None:
            fig, ax = plt.subplots(1,1)
        ax.plot(minmax_hist)
        ax.plot(xs, minmax_dens)
    # use simple neighbor-comparison & peak properties to estimate n peaks
    peaks, props = find_peaks(minmax_dens,
                              height=peak_minheight,
                             )
    npk = len(peaks)
    # calculate absolute lag-1 temporal autocorrelation
    acorr = np.abs(np.corrcoef(hist[:-1], hist[1:])[0,1])
    return npk, acorr, minmax_dens


def estimate_n_peaks(hist,
                     kde_bandwidth=5,
                     peak_minheight=0.6,
                     n_perm=100,
                     plot=False,
                     return_dens=False,
                    ):
    '''
    combine a few metrics to estimate the number of peaks in a histogram
    '''
    # rotate and minmax scale the histogram
    hist, cutidx = rotate_time_series_to_min(hist, return_cutidx=True)
    npk, acorr, scaled_dens = calc_n_kde_peaks(hist,
                                               kde_bandwidth=kde_bandwidth,
                                               peak_minheight=peak_minheight,
                                               plot=False,
                                              )
    # permute the histogram n_perm times, run the same process, collect results
    perm_npk = []
    perm_acorr = []
    for i in range(n_perm):
        hist_perm = np.random.choice(hist, len(hist), replace=False,
                                           )
        assert len(hist_perm) == len(hist)
        pnpk, pacorr, pscaled_dens=calc_n_kde_peaks(hist_perm,
                                                    kde_bandwidth=kde_bandwidth,
                                                    peak_minheight=peak_minheight,
                                                    plot=False,
                                                   )
        perm_npk.append(pnpk)
        perm_acorr.append(pacorr)
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(1,2,1)
        _ = calc_n_kde_peaks(hist,
                             kde_bandwidth=kde_bandwidth,
                             peak_minheight=peak_minheight,
                             plot=True,
                             ax=ax,
                            )
        ax.set_title(f'true data (n peaks={npk})')
        ax = fig.add_subplot(1,2,2)
        ax.hist(perm_acorr, alpha=0.5)
        ax.plot(acorr, 0, 'or')
        ax.set_title('true lag-1 autocorrelation vs. nulls')
    # get empirical P-val from true and permuted absolute lag-1 autocorr vals
    pval = np.mean(np.array(perm_acorr)>= acorr)
    if return_dens:
        scaled_dens_unrotated = np.concatenate((scaled_dens[-cutidx:],
                                                scaled_dens[:-cutidx]))
        return npk, pval, scaled_dens_unrotated
    else:
        return npk, pval



#------------------------------------------------------------------------------
# get histograms and observations
#------------------------------------------------------------------------------
# overall results data struct 
res_dict ={'tid': [],
           'name': [],
           'hist_count': [],
           'hist_count_expec': [],
           'obs_count': [],
           'npeaks': [],
           'npeaks_pval': [],
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
        print((f"\n\n{'.'*80}\nprocessing {tax_name} (TID: {tid}) "
               f"({i+1} of {len(taxa)})..."))

        # get observation histogram
        hist_filename = os.path.join(data_dir,
                        f"hist_data/TID_{tid}_{tax_name.replace(' ', '_')}.csv")
        if not os.path.isfile(hist_filename):
            print("\n\n\thistogram file does not exist; hitting API...\n\n")
            # set var to indicate we need to save the hist fig
            save_hist_fig = True
            hist_success = False
            n_trys=1
            while not hist_success:
                try:
                    # add a 1-second wait, just to help avoid 429 errors
                    time.sleep(1)
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
                    print((f"\n\t\t{np.sum([*hist.values()])} observations returned "
                           f"(vs. {taxa[taxa['tid'] == tid]['count'].values[0]} "
                            "expected based on initial taxa table).\n"))
                except Exception as e:
                    print((f"Error thrown during histogram API call: {e}\n\n"
                           f"Waiting {120 * n_trys} sec, then trying again...\n"))
                    time.sleep(120 * n_trys)
                    n_trys += 1

            # if hist is all zeros then skip and return missing data
            hist_vals = [*hist.values()]
            hist_df = pd.DataFrame({'n_obs': hist_vals}).to_csv(hist_filename,
                                                                index=False)
        else:
            hist_vals = pd.read_csv(hist_filename).loc[:, 'n_obs'].values
            print((f"\n\n\treading histogram of {np.sum(hist_vals)} "
                   "sightings from file...\n"))
            # set var to indicate we don't need to save the hist fig
            if replot_hists:
                save_hist_fig = True
            else:
                save_hist_fig = False
        if np.sum(hist_vals) == 0:
            res_dict['tid'].append(tid)
            res_dict['name'].append(tax_name)
            res_dict['hist_count'].append(np.nan)
            res_dict['hist_count_expec'].append(row['count'])
            res_dict['obs_count'].append(np.nan)
            res_dict['npeaks'].append(np.nan)
            res_dict['npeaks_pval'].append(np.nan)
            res_dict['color'].append(np.nan)
            res_dict['geometry'].append(MultiPolygon())

            # store and save taxon observations
            obs_dict['datetime'].extend([pd.NaT])
            obs_dict['doy'].extend([np.nan])
            obs_dict['doy_circ'].extend([np.nan])
            obs_dict['geometry'].extend([Point()])

        else:
            # get actual observations
            obs_filename = os.path.join(data_dir,
                        f"obs_data/TID_{tid}_{tax_name.replace(' ', '_')}.json")
            if not os.path.isfile(obs_filename):
                print("\tobservation file does not exist; hitting API...\n\n")
                obs_success = False
                n_trys = 1
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
                                print("\t\tgetting observations, page 1...")
                            else:
                                print((f"\t\tgetting observations, page {curr_page} "
                                       f"of {obs_pg_ct}..."))
                            # add a 1-second wait, just to help avoid 429 errors
                            time.sleep(1)
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
                               f"Waiting {120 * n_trys} sec, then trying again...\n"))
                        time.sleep(120 * n_trys)
                        n_trys += 1

                # store taxon observations
                obs_dict['datetime'].extend([str(d).replace(' ', 'T') for d in dates])
                obs_dict['doy'].extend(doys)
                obs_dict['doy_circ'].extend(doys_circ)
                obs_dict['geometry'].extend([Point(c) for c in coords])

                # save the observation data
                obs_df = pd.DataFrame.from_dict(obs_dict)
                obs_gdf = gpd.GeoDataFrame(obs_df,
                                   geometry='geometry',
                                   crs=4326,
                                  )
                obs_gdf.to_file(obs_filename)

            else:
                # read in the observation data instead, if it already exists
                obs_gdf = gpd.read_file(obs_filename)
                print((f"\treading {len(obs_gdf)} observations "
                   "from file...\n\n"))
                coords = [(g.x, g.y) for g in obs_gdf.geometry]

            # estimate the number of peaks in the histogram
            # and the empirical P-value associated with that metric
            npeaks, npeaks_pval, kde = estimate_n_peaks(hist_vals,
                                                  kde_bandwidth=npeak_bandwidth,
                                                  peak_minheight=npeak_minheight,
                                                  n_perm=npeak_nperm,
                                                  plot=False,
                                                  return_dens=True,
                                                  )
            print((f"\n\n\tn fitted peaks: {npeaks} "
                   f"(P-value: {np.round(npeaks_pval, 2)})\n\n"))

            # calculate an alpha hull around the observation coordinates
            hull =  alphashape.alphashape(np.array(coords), alpha=alpha)
            if not isinstance(hull, MultiPolygon):
                hull = MultiPolygon([hull])
            assert isinstance(hull, MultiPolygon)

            # assign a random color to this taxon
            color = get_hex_color()

            # plot hist and both regressions
            fig = plt.figure(figsize=(4,4))
            ax = fig.add_subplot(1,1,1)
            # TODO: FIGURE OUT KDE-UNROTATIING BUG
            ax.plot([*range(len(hist_vals))], hist_vals, ':k')
            ax.plot(np.linspace(0, len(hist_vals), len(kde)),
                    kde*np.max(hist_vals), '-', color=color)
            ax.text(45, 0.9*ax.get_ylim()[1],
                    f"P={np.round(npeaks_pval, 2)}", size=10)
            ax.set_xlabel('week of year')
            ax.set_ylabel('number of observations')
            ax.set_title((f"TID {tid}: {tax_name}\n"
                f"({np.sum(hist_vals)} total flowering observations)"))
            fig_filename = os.path.join(data_dir,
                       f"hist_plots/TID_{tid}_{tax_name.replace(' ', '_')}.png")
            if save_hist_fig:
                print("\tsaving histogram plot...\n")
                fig.savefig(os.path.join(phf.FIGS_DIR,
                                         'inat_phen_hists',
                                         fig_filename),
                            dpi=400)
            plt.close('all')

            # store results
            res_dict['tid'].append(tid)
            res_dict['name'].append(tax_name)
            res_dict['hist_count'].append(np.sum(hist_vals))
            res_dict['hist_count_expec'].append(row['count'])
            res_dict['obs_count'].append(len(coords))
            res_dict['npeaks'].append(npeaks)
            res_dict['npeaks_pval'].append(npeaks_pval)
            res_dict['color'].append(color)
            res_dict['geometry'].append(hull)

        # save what data we have thus far
        res_df = pd.DataFrame.from_dict(res_dict)
        res_gdf = gpd.GeoDataFrame(res_df,
                                   geometry='geometry',
                                   crs=4326,
                                  )
        assert len(np.unique(res_gdf['tid'])) == len(res_gdf)
        res_gdf_curr_it = res_gdf[res_gdf['name'] == tax_name]
        assert len(res_gdf_curr_it) == 1
        assert res_gdf_curr_it['tid'].values[0] == tid
        res_gdf_curr_it.to_file(processed_taxa_filename,
                                mode=processed_file_writemode)
        # switch to appending to the results file after the first taxon,
        # if the script started a new results file and ∴ started in write mode
        if processed_file_writemode == 'w' and i == 0:
            processed_file_writemode = 'a'

        # clear cache every 10 calls, so laptop memory doesn't top out
        if i % 10 == 0:
            print("\n\n\tCLEARING CACHE...\n\n")
            clear_cache()

        stop_time = time.time()
        time_diff = stop_time - start_time
        runtimes.append(time_diff)
        print(("\n\trolling-average runtime: "
               f"{np.round(np.sum(runtimes)/(i+1), 1)} sec per taxon\n"))

    except Exception as e:
        print(f"EXCEPTION THROWN: {e}")
        total_failure = True
        break


#------------------------------------------------------------------------------
# export results
#------------------------------------------------------------------------------
if total_failure:
    print("\n\nTOTAL FAILURE! Interim results saved. Debug and rerun.\n\n")
else:
    if len(taxa) > 0:
        fig, ax = plt.subplots(1)
        res_gdf.to_crs(3857).plot(alpha=0.7,
                                  color=res_gdf['color'],
                                  edgecolor='black',
                                  markersize=10,
                                  ax=ax,
                                  legend='name',
                                 )
        try:
            ctx.add_basemap(ax=ax,
                            source=xyzservices.providers.OpenStreetMap['Mapnik'])
        except Exception as e:
            print("\n\n\tCONTEXTILY ERROR; NO BASEMAP\n\n")
            world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
            world = world[world['continent'] != 'Antarctica']
            world.to_crs(3857).plot(color='none', edgecolor='black', linewidth=0.5,
                                    alpha=0.5, ax=ax, zorder=0)
        fig.savefig(os.path.join(phf.FIGS_DIR, "inat_flower_phen_obs_locs.png"), dpi=600)
    else:
        print("\n\nNo taxa remaining to be processed!\n")

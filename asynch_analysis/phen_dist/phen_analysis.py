import MMRR
import helper_fns as hf
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from datetime import datetime
from shapely.geometry import Point
from greatcircle import GreatCircle
import seaborn as sns
from copy import deepcopy
import re

################
# TODO:

    # BUGS:
        #1. something makes MMRR signif pretty much always...
            # FIX: numpy arrays, when subsetted on a list of row nums and a
            #      list of col nums, pairs them elementwise and extracts
            #      the vector of values at each coord pair,
            #      whereas R takes the rows and cols separately,
            #      giving a new, permuted matrix;
            #      it seems as though this behavior must have changed in numpy,
            #      right?!? becuase I know I tested my MMRR.py script before
            #      and it produced identical results to Ian's MMRR.R...

    # WHY TF DID MMRR RETURN A SIGNIF RESULT (p = 0.02...) ON 2 EQUAL-SIZE
    # MATRICES OF PURELY RANDOM NUMBERS??

    # figure out error involving lat/lon - cell crosswalk code
    # for some of the test spps (it said something about indexing something
    # of length 0, if I remember right?...)

    # check out and scatterplot a handful of test spps

    # further refine the set of spps to be analyzed!
        # after doing that, then decide if the following quesitons are worth
        # attending to:
            # set up to run in parallel on Savio?
            # can I speed up pairwise dist calculations somehow?

    # is it worthwhile to figure out how to try to drop all non-native-range
    # samples?

    # need to reproduce and include the geodetic datum

    # possibly produce a less-filtered map of coefficients,
    # so as to not drop so many samples that fall in no-coefficient cells?

    # will need to manually scour the species, figure out what their native
    # ranges should look like, figure out if they're likely to be domestic,
    # and then drop species and/or samples as necessary

################

# overarching params to control behavior

# set number of species to be analyzed
n_spps_to_analyze = 200
# set limit of number of samples per spp
n_samps_per_spp = 1000
#n_samps_per_spp = None

# set minimum sample size for a valid test
min_samp_size = 30

# whether or not to save the individual scatterplots
savescat=True


# help fns

def scatter_single_spp(phen_dist, seas_dist, geog_dist, spp, savefig=False):
    """
    takes the three pairwise distance matrices for a species
    and returns a scatterplot of them
    """
    # get the indices
    tri_inds = np.tril_indices(phen_dist.shape[0])
    # pull the vals
    phen = phen_dist[tri_inds]
    seas = seas_dist[tri_inds]
    geog = geog_dist[tri_inds]

    # make the plot
    fig, axs = plt.subplots(1,2)
    fig.suptitle('%s:\nphenological dist as a fn of:' % spp, size=20)
    ax1, ax2 = axs
    ax1.set_title('seasonal dist', size=16)
    ax2.set_title('geographic dist', size=16)
    ax1.scatter(seas, phen)
    ax2.scatter(geog, phen)
    ax1.set_ylabel('phen', size=14)
    ax1.set_xlabel('seas', size=14)
    ax2.set_ylabel('phen', size=14)
    ax2.set_xlabel('geog', size=14)
    ax1.set_ylim(0,1)
    ax2.set_ylim(0,1)
    if savefig:
        formatted_spp = re.sub(' ', '_', spp)
        fig.savefig('./phen_dist/results/spp_scat_%s.png' % formatted_spp)
    else:
        return fig


# run analysis
print('\nloading and prepping data\n')

# location of global coefficients raster
coeffs_file = './data/SIF_coeffs.tif'

# load phen data
phen_raw = pd.read_csv('./data/gbif_plant.csv')
#small_phen = pd.read_csv('./data/subset_gbif_data.csv')
#spp_subset = [*small_phen.species.unique()]
#phen_raw = phen_raw[phen_raw.species.isin(spp_subset)]

# only keep the needed columns, then delete raw data to reduce mem usage
phen = phen_raw.loc[:, ['species',
                        'decimalLongitude',
                        'decimalLatitude',
                        'eventDate']]
del phen_raw

# get number of species in dataset
n_spps = len(phen.species.unique())

# drop nas
for col in ['eventDate', 'species', 'decimalLatitude', 'decimalLongitude']:
    phen = phen[phen[col].notna()]

# convert day to doy in radians
# TODO: decide if leap years might create a problem
rad_day = [datetime.fromisoformat(row['eventDate']).timetuple(
            ).tm_yday/365*2*np.pi for n, row in phen.iterrows()]
phen.loc[:, 'rad_day'] = rad_day

# dict to hold MMRR output
results = {col:[] for col in ['spp', 'n', 'R^2', 'F-statistic', 'F p-value',
                              'Intercept', 'Intercept(t)', 'Intercept(p)',
                              'seas', 'seas(t)', 'seas(p)',
                              'geog', 'geog(t)', 'geog(p)']}

# dict to hold all the figures
figs = {}

# loop over species
spp_n = 0
for spp in phen.species.unique():
    if spp_n < n_spps_to_analyze:
        try:
            print(('\n now analyzing species %i of %i to be analyzed '
                   '(%i total spps):'
                   '\n\t%s\n') % (spp_n, n_spps_to_analyze, n_spps, spp))

            # subset for focal species
            spp_phen = phen[phen.species == spp]
            print('\n\thas %i samples\n' % len(spp_phen))

            # keep only rows with unique site localities,
            # getting the mean date (in radians) for each locality
            # (to avoid division-by-zero problems in MMRR, because identical
            #  sites give identical seasonal and geographic distance columns,
            #  creating perfect multicollinearity in the OLS models)
            # NOTE: resetting index so lat and lon to behave like cols again
            spp_phen = spp_phen.loc[:, ['decimalLatitude',
                                        'decimalLongitude',
                                        'rad_day']].groupby(['decimalLongitude',
                                    'decimalLatitude']).mean().reset_index()
            print(('\n\thas %i samples after '
                   'dropping duplicate sites\n') % len(spp_phen))

            # coerce to a GeoDataFrame
            # TODO: need to pull the data again
            #       and include the geodeticDatum column!
            #       for now, assuming EPSG:4326,
            #       per the note on GBIF's website for
            #       any rows with no data in that column anyhow;
            #       https://www.gbif.org/data-quality-requirements-occurrences
            pts = [Point(row.decimalLongitude,
                         row.decimalLatitude) for n, row in spp_phen.iterrows()]
            spp_phen = gpd.GeoDataFrame(spp_phen, geometry=pts, crs='EPSG:4326')

            # reset the index (to avoid row-number problems)
            spp_phen.reset_index(inplace=True)

            # cut down the size of the dataframe, if required
            if n_samps_per_spp is not None and len(spp_phen) > n_samps_per_spp:
                keep_rows = np.random.choice([*range(len(spp_phen))],
                                             n_samps_per_spp, replace=False)
                spp_phen = spp_phen.iloc[keep_rows, :].reset_index()
                print(('\n\thas %i samples after '
                   'limiting number of samples\n') % len(spp_phen))

            # calculate seas dist matrix
            print('\n\t calculating seas_dist matrix\n')
            pt_coords = spp_phen.loc[:, ['decimalLongitude',
                                         'decimalLatitude']].values
            seas_dist = hf.get_seasonality_info_points(coeffs_file, pt_coords)
            # will drop rows and cols for any samples
            # that fall in cells without fitted coeffs
            not_missing_seas_dist = np.where(np.nansum(seas_dist, axis=0)>0)[0]

            # create dist matrices to be filled up
            phen_dist = np.zeros([len(spp_phen)]*2) * np.nan
            geog_dist = np.zeros([len(spp_phen)]*2) * np.nan

            # loop over sample pairs and calculate their distances
            print('\n\t calculating phen_dist and geog_dist matrices\n')
            for i, row in spp_phen.iterrows():
                if i in not_missing_seas_dist:
                    # set diag to 0
                    phen_dist[i,i] = 0
                    geog_dist[i,i] = 0
                    for j in range(i+1, len(spp_phen)):
                        if j in not_missing_seas_dist:

                            # set dist to 0 if coords are identical
                            if (row.decimalLatitude == spp_phen.iloc[j,
                                                :]['decimalLatitude'] and
                                row.decimalLongitude == spp_phen.iloc[j,
                                                :]['decimalLongitude']):
                                dist = 0
                            else:
                                # calculate the geog dist (great circle dist)
                                gc = GreatCircle()
                                gc.latitude1_degrees = row.decimalLatitude
                                gc.longitude1_degrees = row.decimalLongitude
                                gc.latitude2_degrees = spp_phen.iloc[j,
                                                        :]['decimalLatitude']
                                gc.longitude2_degrees = spp_phen.iloc[j,
                                                        :]['decimalLongitude']
                                gc.calculate()
                                dist = gc.distance_kilometres*1000
                            # add to geog_dist matrix
                            geog_dist[i,j] = dist
                            geog_dist[j,i] = dist

                            # calculate the phenological distance
                            # (chord dist on unit circle)
                            ang_diff = row.rad_day - spp_phen.iloc[j,
                                                                   :]['rad_day']
                            dist = np.sin(np.abs(ang_diff)/2)
                            # add to geog_dist matrix
                            phen_dist[i,j] = dist
                            phen_dist[j,i] = dist

            # calculate seas dist matrix
            #print('\n\t calculating seas_dist matrix\n')
            #pt_coords = spp_phen.loc[:, ['decimalLongitude',
            #                             'decimalLatitude']].values
            #seas_dist = hf.get_seasonality_info_points(coeffs_file, pt_coords)

            #assert(phen_dist.shape == geog_dist.shape == seas_dist.shape), (''
            #                                'DIST MATRIX SHAPES NOT EQUAL!')

            # drop rows and cols for any samples
            # that fall in cells without fitted coeffs
            not_missing = np.where(np.nansum(seas_dist, axis=0)>0)[0]
            phen_dist = phen_dist[:, not_missing][not_missing,:]
            seas_dist = seas_dist[:, not_missing][not_missing,:]
            geog_dist = geog_dist[:, not_missing][not_missing,:]

            print(('\n\thas %i samples after dropping '
                   'those without seasonal distance '
                   'values\n') % seas_dist.shape[0])

            assert(phen_dist.shape == geog_dist.shape == seas_dist.shape), (''
                                            'DIST MATRIX SHAPES NOT EQUAL!')


            if len(not_missing) < min_samp_size:
                raise Exception('sample size of %i inadequate.' % len(
                                                                not_missing))

            # run MMRR
            mod = MMRR.MMRR(phen_dist, [seas_dist, geog_dist], ['seas', 'geog'])

            # save the results
            results['spp'].append(spp)
            results['n'].append(seas_dist.shape[0])
            for k,v in mod.items():
                results[k].append(v)

            # produce the scatterplot
            figs[spp] = scatter_single_spp(phen_dist, seas_dist, geog_dist,
                                           spp, savefig=savescat)
            spp_n += 1

        # if the species threw an error, report it, then move on
        except Exception as e:
            print('\nERROR FOR SPECIES %s:\n\t%s\n' % (spp, e))
            results['spp'].append(spp)
            results['n'].append(np.nan)
            for k in results.keys():
                if k not in ['spp', 'n']:
                    results[k].append(np.nan)


# cast results as DataFrame and save
results_df = pd.DataFrame.from_dict(results)
results_df.to_csv('phen_results.csv', index=False)

# violin plot of results
fig, axs = plt.subplots(1,2)
ax1, ax2 = axs
# plot p-values
viol_df = results_df.melt('spp', value_vars=['geog(p)', 'seas(p)', 'F p-value'])
#viol_df['value'] = np.log(viol_df['value'])
sns.violinplot(x='variable', y='value', data=viol_df, ax=ax1,
               inner=None, palette='Set3')
ax1.set_ylabel('p-value', size=18)
#ax1.set_ylabel('log p-value')
ax1.set_xlabel('stat', size=18)
ax1.tick_params(labelsize=14)
ax1.plot([-1,5], [0.05]*2, ':k')
ax1.set_xlim(-1,3)
# plot coeff values
viol_df = results_df.melt('spp', value_vars=['geog', 'seas'])
sns.violinplot(x='variable', y='value', data=viol_df, ax=ax2,
               inner=None, palette='Set3')
ax2.set_ylabel('coefficient value', size=18)
ax2.set_xlabel('coefficient', size=18)
ax2.tick_params(labelsize=14)
ax2.plot([-1,6], [0]*2, ':k')
ax2.set_xlim(-1,2)
fig.show()


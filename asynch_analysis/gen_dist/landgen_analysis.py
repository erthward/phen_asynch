import helper_fns as hf
from MMRR import MMRR

import os
import glob
import numpy as np
import pandas as pd
import seaborn as sns
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from skbio.stats.distance import mantel


#############################################################################
# TODO:

    # fix map code

    # what to do about having so few samples for pop-based analysis?
    #   could probably write some code to draw seasonality from nearest 
    #   non-missing cell, within some tolerance distance?

    # is it a problem if permuted matrix equal to original in small-n MMRR?
    # maybe just toss the matrix and try again, rather than throwing assertion
    # error?

    # decide whether log-transform of matrices makes sense (a la Kling);
    # perhaps Shapiro tests could help determine how frequently the data looks
    # markedly non-normal?
    # but first and foremost, I don't believe normality is an assumption of a
    # Mantel test anyhow?...

    # 2x-check/figure out Smouse genetic dist metric (use other?) (how not drop
    # some pops?)

    # correct for multiple tests

    # mean coeff in sig vs insig seasonality spps

    # mean asynch within convex hull of points for significant-seasonality spps
    #     vs seasonality-inisigif spps?

    # map pops colored by: both insig, clim insig, seas insig, both sig
#############################################################################


# MAIN ANALYSIS PARAMETERS

# genetic distance metric to use
#gen_dist_metric = 'smouse'
gen_dist_metric = 'wc84'

# max horizontal and vertical 'radius' of a nan-patch in the asynch raster
# that will be filled by interpolation using rasterio.fill.fillnodata
fill_tol = 5

# decisions on how to handle the distance-matrix data
standardize_dist_mats = True
log_transform_dist_mats = False


# DATA/FILE PARAMS
gen_subdir = 'kling_gendata'
genepop_file_patt = "%s_genepop.txt"
gen_dist_file_patt = "%s_genepop.gen." + gen_dist_metric + ".csv"
csv_file_patt = "%s_populations.csv"


# HELPER FNS
def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


# DATA PREP

# get a list of all the species
spps = [os.path.basename(f).split('_')[0] for f in glob.glob(
                    os.path.join(hf.DATA_DIR, gen_subdir + '/*_genepop.gen'))]
spps = [*np.unique(spps)]

# data structures to hold MMRR results and sample sizes
res = {}
ns = {}
mantel_coeffs = {}
mantel_ps = {}
sea_geo_corrs = {}
sea_env_corrs = {}


# SPECIES ANALYSES
for spp in spps:

    # try running the landgen analysis, if the genetic distance matrix
    # file exists for this species
    exists = os.path.isfile(os.path.join(hf.DATA_DIR, gen_subdir,
                                         gen_dist_file_patt % spp))
    if not exists:
        print('\nNO GENETIC DISTANCE MATRIX FILE FOR SPECIES %s\n' % spp)
    else:
        print('\nNOW PROCESSING SPECIES %s\n' % spp)

        try:

            # get gen distmat
            gen_df = pd.read_csv(os.path.join(hf.DATA_DIR, gen_subdir,
                                              gen_dist_file_patt % spp))
            gen = gen_df.iloc[:, 1:].values

            # get pops' coords
            pops = pd.read_csv(os.path.join(hf.DATA_DIR, gen_subdir,
                                                  csv_file_patt % spp))

            # get indivds' coords from pops' and individs' pops
            pops['pop'] = ['p%i' % i for i in pops['ID']]
            individs_pops = [s.split('_')[0] for s in gen_df.iloc[:,0]]
            merged = pd.merge(pd.DataFrame({'pop':individs_pops}),
                              pops, on='pop')
            pts= merged.loc[:, ['longitude', 'latitude']].values

            # get geo distmat
            geo = hf.calc_pw_geo_dist_mat(pts)

            # get environmental (i.e. bioclim) distmat
            env = hf.calc_pw_clim_dist_mat(pts)

            # get seasonal distmat
            #sea = hf.get_raster_info_points(hf.COEFFS_FILE, pts, 'ts_pdm',
            sea = hf.get_raster_info_points(hf.COEFFS_FILLED_FILE, pts, 'ts_pdm',
                                            fill_nans=False, fill_tol=fill_tol)

            # drop pts for which there are no seasonal distance values, no
            # environmental distance values, or no genetic distance values
            not_missing = np.where(np.nansum(sea, axis=0)>0)[0]
            sea = sea[:, not_missing][not_missing,:]
            env = env[:, not_missing][not_missing,:]
            geo = geo[:, not_missing][not_missing,:]
            gen = gen[:, not_missing][not_missing,:]
            again_not_missing = np.where(np.nansum(env, axis=0)>0)[0]
            sea = sea[:, again_not_missing][again_not_missing,:]
            env = env[:, again_not_missing][again_not_missing,:]
            geo = geo[:, again_not_missing][again_not_missing,:]
            gen = gen[:, again_not_missing][again_not_missing,:]
            still_not_missing = np.where(np.nansum(gen, axis=0)>0)[0]
            sea = sea[:, still_not_missing][still_not_missing,:]
            env = env[:, still_not_missing][still_not_missing,:]
            geo = geo[:, still_not_missing][still_not_missing,:]
            gen = gen[:, still_not_missing][still_not_missing,:]

            # check that matrices are symmetric and hollow
            for matname, mat in {'gen': gen,
                                 'sea': sea,
                                 'env': env,
                                 'geo': geo}.items():
                assert np.all(mat[np.diag_indices(mat.shape[0])] == 0), ('mat'
                                            'rix %s not hollow!') % mat_name
                assert check_symmetric(mat), ('matrix %s '
                                              'not symmetric!') % mat_name

            # get sample size
            n = gen.shape[0]
            print('\n\tSpecies %s has %i samples\n' % (spp, n))

            # get covariance between geo and sea
            sea_geo_corr = np.corrcoef(geo[np.triu_indices(geo.shape[0])],
                                       sea[np.triu_indices(sea.shape[0])])[1,0]
            sea_env_corr = np.corrcoef(env[np.triu_indices(env.shape[0])],
                                       sea[np.triu_indices(sea.shape[0])])[1,0]
            sea_geo_corrs[spp] = sea_geo_corr
            sea_env_corrs[spp] = sea_env_corr

            # run mantel between gen and sea
            mantel_coeff, mantel_p, mantel_n = mantel(gen, sea)
            mantel_coeffs[spp] = mantel_coeff
            mantel_ps[spp] = mantel_p

            if log_transform_dist_mats:
                gen = np.log(gen)
                sea = np.log(sea)
                env = np.log(env)
                geo = np.log(geo)

            if standardize_dist_mats:
                # standardize all the arrays
                gen = hf.standardize_array(gen)
                geo = hf.standardize_array(geo)
                env = hf.standardize_array(env)
                sea = hf.standardize_array(sea)


            # run MMRR
            mod = MMRR(gen, [sea, env, geo], ['sea', 'env', 'geo'])

            # store the results
            res[spp] = mod
            ns[spp] = n

        except Exception as e:
            print('\nERROR THROWN FOR SPECIES %s\n\t%s\n' % (spp, e))
            import traceback
            traceback.print_exc()

# digest all results into a DataFrame
merged_res_dict = {k:[] for k in [*res.values()][0].keys()}
merged_res_dict['spp'] = []
merged_res_dict['n'] = []
merged_res_dict['mantel_coeff'] = []
merged_res_dict['mantel_p'] = []
merged_res_dict['sea_geo_corr'] = []
merged_res_dict['sea_env_corr'] = []
for k, res_subdict in res.items():
    merged_res_dict['spp'].append(k)
    merged_res_dict['n'].append(ns[k])
    merged_res_dict['mantel_coeff'].append(mantel_coeffs[k])
    merged_res_dict['mantel_p'].append(mantel_ps[k])
    merged_res_dict['sea_geo_corr'].append(sea_geo_corrs[k])
    merged_res_dict['sea_env_corr'].append(sea_env_corrs[k])
    for sub_k, sub_v in res_subdict.items():
        merged_res_dict[sub_k].append(sub_v)
results_df = pd.DataFrame.from_dict(merged_res_dict)
results_df.to_csv('landgen_results.csv', index=False)

# violin plot of results
inner='box'
fig, axs = plt.subplots(1,2)
ax1, ax2 = axs
# plot p-values
viol_df = results_df.melt('spp', value_vars=['geo(p)', 'env(p)', 'sea(p)',
                                             'F p-value'])
sns.violinplot(x='variable', y='value', data=viol_df, ax=ax1,
               inner=inner, palette='Set3')
ax1.set_ylabel('p-value', size=18)
#ax1.set_ylabel('log p-value')
ax1.set_xlabel('stat', size=18)
ax1.tick_params(labelsize=14)
ax1.plot([-1,5], [0.05]*2, ':k')
ax1.set_xlim(-1,4)
# plot coeff values
viol_df = results_df.melt('spp', value_vars=['geo', 'env', 'sea'])
sns.violinplot(x='variable', y='value', data=viol_df, ax=ax2,
               inner=inner, palette='Set3')
ax2.set_ylabel('coefficient value', size=18)
ax2.set_xlabel('coefficient', size=18)
ax2.tick_params(labelsize=14)
ax2.plot([-1,6], [0]*2, ':k')
ax2.set_xlim(-1,4)
fig.show()


# create colors for the 4 types of model significance
# NOTE: keys are structured as (env sig?, sea sig?)
both_sig =    "#c542f5" # purple = both sig
sea_sig =     "#f54257" # red = sea sig
env_sig =     "#42b3f5" # blue = env sig
both_insig = "#000000" # black = both insig
mod_insig =   "#bdbdbd" # grey = model itself insig

sig_cols = {(True, True):   both_sig,
            (False, True):  sea_sig,
            (True, False):  env_sig,
            (False, False): both_insig,
            None:           mod_insig,
           }

leg_cols = {'both':      both_sig,
            'sea':       sea_sig,
            'env':       env_sig,
            'neither':   both_insig,
            'mod insig': mod_insig,
           }

# plot the countries map
countries = gpd.read_file(os.path.join(hf.COUNTRIES_DATA_DIR, 'countries.shp'))

fig_map, ax = plt.subplots(1,1)
countries.plot(ax=ax, color='white', edgecolor='black', linewidth=0.3)

for spp in spps:

    print('\nNOW MAPPING SPECIES %s\n' % spp)

    if spp in [*results_df['spp']]:
        # get pops' coords
        pops = pd.read_csv(os.path.join(hf.DATA_DIR, gen_subdir,
                                        csv_file_patt % spp))

        results_subdf = results_df[results_df['spp'] == spp]
        assert len(results_subdf) == 1, ('more than one results row for '
                                         'species %s!') % spp
        if results_subdf['F p-value'].values[0] >= 0.05:
            col_code = None
        else:
            col_code = (results_subdf['env(p)'].values[0] < 0.05,
                        results_subdf['sea(p)'].values[0] < 0.05)

        # get spp's color
        col = sig_cols[col_code]

        # scatter the pops
        ax.scatter(x=pops['longitude'].iloc[0], y=pops['latitude'].iloc[0],
                   #s=results_subdf['n'].iloc[0],
                   c=col, edgecolor='k', linewidth=0.5)

# add the legend
legend_elements = [Line2D([0], [0],
                          marker='o', color= 'w',
                          label=l, markerfacecolor=c,
                          markeredgecolor='k', linewidth=0.5,
                          markersize=7) for l, c in leg_cols.items()]
ax.legend(handles=legend_elements, loc='lower left')
fig_map.show()


# plot boxplots split by whether or not a mantel test of gen ~ sea is sig
fig_box, ax_box = plt.subplots(1,1)
results_df['mantel_p_lt_p05'] = results_df['mantel_p'] < 0.05
results_df.boxplot(column=['env(p)', 'geo(p)', 'sea(p)',
                           'sea_env_corr', 'sea_geo_corr', 'mantel_coeff',
                           'F p-value', ],
                   by = 'mantel_p_lt_p05', ax=ax_box)
fig_box.show()



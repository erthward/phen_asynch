import os
import numpy as np
import pandas as pd
import pyinaturalist as pynat
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, find_peaks_cwt, peak_prominences, peak_widths
from sklearn.neighbors import KernelDensity
from sklearn.mixture import GaussianMixture


# TODO:

    # consider just calculating temporal autocorrelation!?

    # finalize and test out on some new taxa

    # integrate into map_....py

    # update map...py to only download obserations if not already saved
    # and to download and saves hists

    # also update map...py to save npeaks and its performance metric(s)

    # rerun!

    # add npeaks/etc into mapping stuf

    # circle back to Ian




# download histograms that I see as clear examples of 0, 1, and 2 peaks
spp = ['Metrosideros polymorpha',  # 0
       'Trillium grandiflorum',    # 1
       'Stephanomeria pauciflora', # 2
      ]
taxa = pd.read_csv('./all_inat_plant_phen_taxa_w_TRY_pgf.csv')
spp = {taxa[taxa['name']==sp]['tid'].iloc[0]: sp for sp in spp}


# max positional accuracy value allowed
# NOTE: accuracy values are expressed as "accurate to within <VAL> meters"
#       (from iNat site: "distance in meters that includes the entire area
#        where the observation could have taken place")
max_pos_acc_val = 1000

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


hists = {}
for tid, sp in spp.items():
    hist = pynat.get_observation_histogram(taxon_id=tid,
                                           term_id=term_id,
                                           term_value_id=term_value_id,
                                           date_file=date_file,
                                           interval=interval,
                                           quality_grade=quality_grade,
                                           native=native,
                                           captive=captive,
                                          )
    hist_vals = [*hist.values()]
    hists[tid] = hist_vals

fig = plt.figure()
ct = 1
for tid, hist in hists.items():


    ax = fig.add_subplot(1,3,ct)
    ax.plot(ts)
    ct+=1
fig.show()



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
    xs = np.arange(0, len(x)+0.1, 0.1).reshape((-1, 1))
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
    # fit a Gaussian mixture model using that number of peaks, save AIC and BIC
    D = minmax_dens.reshape((-1, 1))
    gmm = GaussianMixture(npk).fit(D)
    aic = gmm.aic(D)
    bic = gmm.bic(D)
    # calculate absolute lag-1 temporal autocorrelation
    acorr = np.abs(np.corrcoef(hist[:-1], hist[1:])[0,1])
    return npk, aic, bic, acorr


def estimate_n_peaks(hist,
                     kde_bandwidth=5,
                     peak_minheight=0.6,
                     n_perm=100,
                     plot=False,
                    ):
    '''
    combine a few metrics to estimate the number of peaks in a histogram
    '''
    # rotate and minmax scale the histogram
    hist = rotate_time_series_to_min(hist)
    npk, aic, bic, acorr = calc_n_kde_peaks(hist,
                                          kde_bandwidth=kde_bandwidth,
                                          peak_minheight=peak_minheight,
                                          plot=False,
                                         )
    # permute the histogram n_perm times, run the same process, and collect AICs
    perm_npk = []
    perm_aic = []
    perm_bic = []
    perm_acorr = []
    for i in range(n_perm):
        hist_perm = np.random.choice(hist, len(hist), replace=False,
                                           )
        assert len(hist_perm) == len(hist)
        pnpk, paic, pbic, pacorr = calc_n_kde_peaks(hist_perm,
                                                 kde_bandwidth=kde_bandwidth,
                                                 peak_minheight=peak_minheight,
                                                 plot=False,
                                                )
        perm_npk.append(pnpk)
        perm_aic.append(paic)
        perm_bic.append(pbic)
        perm_acorr.append(pacorr)
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(2,2,1)
        _ = calc_n_kde_peaks(hist,
                             kde_bandwidth=kde_bandwidth,
                             peak_minheight=peak_minheight,
                             plot=True,
                             ax=ax,
                            )
        ax.set_title(f'true data (n peaks={npk})')
        ax = fig.add_subplot(2,2,2)
        ax.hist(perm_aic, alpha=0.5)
        ax.plot(aic, 0, 'or')
        ax.set_title('true AIC vs. nulls')
        ax = fig.add_subplot(2,2,3)
        ax.hist(perm_bic, alpha=0.5)
        ax.plot(bic, 0, 'or')
        ax.set_title('true BIC vs. nulls')
        ax = fig.add_subplot(2,2,4)
        ax.hist(perm_acorr, alpha=0.5)
        ax.plot(acorr, 0, 'or')
        ax.set_title('true lag-1 autocorrelation vs. nulls')
    # get empirical P-val from true and permuted absolute lag-1 autocorr vals
    pval = np.mean(np.array(perm_acorr)>= acorr)
    return npk, pval



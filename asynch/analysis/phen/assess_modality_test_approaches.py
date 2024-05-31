import pandas as pd
import numpy as np
import diptest
import matplotlib.pyplot as plt
from scipy.stats import wasserstein_distance

#----------------------------------------------------------------------
# dip test:
#----------------------------------------------------------------------

# read in examples of time series with the following numbers of obvious
# modes in their iNaturalist flowering phenology histograms:
#     0 TID: 99999: this is just simulated data, ~N(17, 6), rounded to ints,
#                   then clipped to minimum 0, since I have yet to notice any
#                   species with no clear seasonal peaks...
#     1 TID: 55402: Trillium grandiflorum
#     2 TID: 58040: Stephanomeria pauciflora

data = pd.read_csv('./example_flow_phen_hists_w_0_1_2_modes.csv')

minmax_scale = lambda a: (a - np.min(a))/(np.max(a) - np.min(a))

fig = plt.figure(figsize=(12,6))
fig.suptitle('dip test', size=15)
colors = ['black', 'blue', 'red']
for i in range(3):
    x = data.iloc[:, i]
    dip, pval = diptest.diptest(x)
    #dipscale, pvalscale = diptest.diptest(minmax_scale(x))
    ax = fig.add_subplot(3,1,i+1)
    ax.plot([*range(len(x))], x, color=colors[i], label=f"{i} modes")
    ax.text(5, 0.3*np.max(x), f"dip: {np.round(dip, 2)} (P: {np.round(pval, 2)})")
    #ax.text(5, 0.15*np.max(x), f"scaled dip: {np.round(dip, 2)} (scaled P: {np.round(pval, 2)})")
fig.legend()
fig.show()
fig.savefig('modality_determination_by_diptest.png', dpi=400)


#-----------------------------------------------------------------------
# earth-mover's distance between time-shifted time series and sin curves
#-----------------------------------------------------------------------

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


# construct canonical 1x and 2x sinusoidal curves, against which data will be
# compared using earth-mover's distance
sin_1x = minmax_scale(np.cos(np.linspace(np.pi, 3*np.pi, 53)))
sin_2x = minmax_scale(np.cos(np.linspace(np.pi, (3+2)*np.pi, 53)))
sin_0x = np.ones((53)) * np.mean([np.mean(s) for s in [sin_1x, sin_2x]])

fig = plt.figure(figsize=(12,6))
colors = ['black', 'blue', 'red']
for i in range(3):
    x = data.iloc[:, i]
    xscale = minmax_scale(x)
    xrot = rotate_time_series_to_min(xscale)
    eud0 = np.linalg.norm(sin_0x-xrot)
    eud1 = np.linalg.norm(sin_1x-xrot)
    eud2 = np.linalg.norm(sin_2x-xrot)
    euds = [eud0, eud1, eud2]
    emd0 = wasserstein_distance(sin_0x, xrot)
    emd1 = wasserstein_distance(sin_1x, xrot)
    emd2 = wasserstein_distance(sin_2x, xrot)
    emds = [emd0, emd1, emd2]
    min_eud = np.argmin(euds)
    min_emd = np.argmin(emds)
    dists = [emds, euds]
    mins = [min_emd, min_eud]
    linetypes = {False: ':', True: '-'}
    for j in range(2):
        ax = fig.add_subplot(3,2,(2*i)+j+1)
        for ns, s in enumerate([sin_0x, sin_1x, sin_2x]):
            ax.plot([*range(len(s))],
                    s,
                    linetypes[ns==mins[j]],
                    color=colors[ns],
                   )
            ax.text(5, (0.4 -0.1*(ns+1))*np.max(xrot),
                        f"dist: {np.round(dists[j][ns], 2)}",
                        color=colors[ns],
                       )
            if ns == mins[j]:
                ax.plot([*range(len(xrot))], xrot, '-', color=colors[ns])
        if i == 0:
            if j == 0:
                ax.set_title('earth-mover\'s dist')
            else:
                ax.set_title('Euclidean dist')
fig.show()
fig.savefig('modality_determination_by_dist.png', dpi=400)



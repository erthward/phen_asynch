import numpy as np
import pandas as pd
import subprocess
import os

def run_diptest_in_R(vals, is_hist=True, noise_sigma=0):
    '''
    use the histogram values to run the dip test in R
    If input vals are samples from the distribution to be tested then set
    is_hist to False.
    If input vals are densities within the subsequent bins of a histogram
    calculated from samples from that distribution then set is_hist to True.
    '''
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


def get_samp(bi=False, n=1000):
    assert n%2 == 0, 'n must be even!'
    mu = np.random.normal(15, 1)
    sigma = np.random.normal(2, 0.5)
    if bi:
        mu = [mu, np.random.normal(30, 1)]
        sigma = [sigma, np.random.normal(2, 0.5)]
        samps = [np.random.normal(mu[i], sigma[i], int(n/2)) for i in range(2)]
        samps = np.array([*samps[0]] + [*samps[1]])
    else:
        samps = np.random.normal(mu, sigma, n)
    assert len(samps) == n
    return samps


# number of simulated samples to test
n_samps = 100
n = 1000
assert n_samps%2 == 0, 'n must be even!'
print((f"\nDrawing and testing {n} samples each of {n_samps} distributions, "
       "a randomly chosen half of which will be bimodal."))

# randomly choose half of the simulated samples to be bimodal
bimod_idxs = np.random.choice([*range(n_samps)], int(n_samps/2), replace=False)

# draw all random samples from either uni- or bimodal normal dists
samps = [get_samp(i in bimod_idxs, n=n) for i in range(n_samps)]

dips = []
ps = []
# run the dip test for each sample and record the dip stats and p-values
for i, samp in enumerate(samps):
    dip, p = run_diptest_in_R(samp, is_hist=False, noise_sigma=0)
    dips.append(dip)
    ps.append(p)

# check sensitivity and specificity
bimod_inferred = np.argwhere(np.array(ps)<0.05).ravel()
sensitivity = np.mean([i in bimod_inferred for i in bimod_idxs])
specificity = (1-(np.mean([i not in bimod_idxs for i in bimod_inferred])))
print('\n\n')
print(f"{np.round(100 * sensitivity, 2)}% sensitivity")
print(f"{np.round(100 * specificity, 2)}% specificity\n\n")

print(f"(true biomdals: {sorted([*bimod_idxs])})")
print(f"(inferred bimodals: {[*bimod_inferred]})\n\n")

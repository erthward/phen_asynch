import numpy as np
import matplotlib.pyplot as plt

t = np.arange(0, 2*np.pi, 0.01)

def reg(t, beta_sin1, beta_cos1, beta_sin2, beta_cos2, var):
    # add 2 pi to t
    t2pi = t + 2*np.pi
    # calculate phase and amplitude vals for each freq
    phase1 = np.arctan2(beta_cos1, beta_sin1)
    amp1 = np.sqrt(beta_sin1**2 + beta_cos1**2)
    phase2 = np.arctan2(beta_cos2, beta_sin2)
    amp2 = np.sqrt(beta_sin2**2 + beta_cos2**2)
    phase_tot = np.sqrt(amp1*(phase1**2)/(amp1+amp2) + amp2*(phase2**2)/(amp1+amp2))
    # calculate true line, and data
    true = (beta_sin1 * np.sin(t2pi) +
           beta_cos1 * np.cos(t2pi) + 
           beta_sin2 * np.sin(2*t2pi) +
           beta_cos2 * np.cos(2*t2pi))
    data = true + np.random.normal(0, var, len(t))
    # get y limits, create plot
    ylim = np.max(np.abs(data))
    fig = plt.figure()
    plt.plot(t, [0]*len(t), ':k', lw=0.2)
    # plot data scattered around true line
    plt.plot(t, true, c='plum')
    plt.scatter(t, data, c='gray', s=0.5)
    # plot freq1 components and phase, if necessary
    if np.abs(beta_sin1) > 0 or np.abs(beta_cos1) > 0:
        plt.plot(t, np.sin(t), '-b', lw=beta_sin1)
        plt.plot(t, np.cos(t), '-b', lw=beta_cos1)
        plt.plot([phase1]*2, [-ylim, ylim], '--b', lw=0.5)
    # plot freq2 components and phase, if necessary
    if np.abs(beta_sin2) > 0 or np.abs(beta_cos2) > 0:
        plt.plot(t, np.sin(2 * t), '-r', lw=beta_sin2)
        plt.plot(t, np.cos(2 * t), '-r', lw=beta_cos2)
        plt.plot([phase2]*2, [-ylim, ylim], '--r', lw=0.5)
    #plot weighted overall phase
    plt.plot([phase_tot]*2, [-ylim, ylim], '--', c='purple', lw=2)
    plt.plot([t for n, t in enumerate(t) if true[n] == max(true)]*2,
             [-ylim, ylim], c='purple', lw=2)
    # set plot limits
    plt.xlim(0, 2*np.pi)
    plt.ylim(-ylim, ylim)
    print('phase 1', phase1, 'phase 2', phase2, 'amp 1', amp1, 'amp 2', amp2,
          'OVERALL PHASE', phase_tot)
    #return (true, data)




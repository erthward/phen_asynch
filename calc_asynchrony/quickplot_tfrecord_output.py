import numpy as np
import matplotlib.pyplot as plt
import asynch_fns as af

def plot_n(n):
    ip, ipnull = af.get_inpatches_outpatches(
        './global_coeffs_SIF-OUT-%s.tfrecord' % str(n).zfill(5),
        af.OUTBANDS, [d-60 for d in af.DIMS]);
    ps = [p for p in ip]
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.imshow(ps[0][0,:,:], interpolation='nearest', cmap='magma')
    ax2 = fig.add_subplot(122)
    ax2.imshow(ps[0][0,:,:], interpolation='nearest', cmap='magma')
    plt.show()



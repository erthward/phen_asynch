import numpy as np
import matplotlib.pyplot as plt
import rioxarray as rxr
import os, re

data_dir = '/media/deth/SLAB/diss/3-phn/corr_data'

for var in ['NIRv', 'SIF']:
    for neigh_rad in [50, 100, 150]:
        print('\n\nPLOTTING ALL SHAP MAPS FOR VAR %s, NEIGH RAD %i KM...\n\n' %
              (var, neigh_rad))
        files = os.listdir(data_dir)
        files = [f for f in files if re.search('^SHAP_map_.*_%s_%ikm.tif' %
                                               (var, neigh_rad), f)]

        fig = plt.figure(figsize=(16,20))
        nrow = 6
        ncol = 2
        for i, f in enumerate(files):
            ax = fig.add_subplot(nrow, ncol, i+1)
            rast = rxr.open_rasterio(os.path.join(data_dir, f), masked=True)
            minmax_val = max(np.abs([np.nanpercentile(rast[0], 1),
                                     np.nanpercentile(rast[0],99)]))
            rast[0].plot.imshow(ax=ax,
                               cmap='coolwarm',
                               vmin=-minmax_val,
                               vmax=minmax_val)
            ax.set_title(os.path.split(f)[-1])
        fig.savefig(os.path.join(data_dir,
                                 'SHAP_maps_%s_%ikm.png' % (var, neigh_rad)),
                                 dpi=700)

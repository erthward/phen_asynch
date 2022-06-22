import matplotlib.pyplot as plt
import geopandas as gpd
import rioxarray as rxr
import numpy as np
import os, re

countries = gpd.read_file(('/home/deth/Desktop/CAL/research/projects/'
                           'seasonality/results/maps/NewWorldFile_2020.shp'))

data_dir = '/media/deth/SLAB/seasonality/results/'
files = [f for f in os.listdir(data_dir) if re.search('_SHAP_vals', f)]

for file in files:
    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(1,1,1)
    rast = rxr.open_rasterio(os.path.join(data_dir, file), masked=True)[0]
    rast = rast.rio.write_crs(4326).rio.reproject(8857)
    rast.plot.imshow(ax=ax,
                     zorder=0,
                     cmap='coolwarm_r',
                     vmin=np.nanpercentile(rast, 0.25),
                     vmax=np.nanpercentile(rast, 99.75),
                     add_colorbar=True,
                    )
    countries.to_crs(8857).plot(ax=ax,
                                color='none',
                                edgecolor='black',
                                linewidth=0.5,
                                alpha=0.8,
                                zorder=1,
                               )
    ax.set_xlim(rast.rio.bounds()[0::2])
    ax.set_ylim(rast.rio.bounds()[1::2])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_title(os.path.splitext(file)[0].split('_')[-1],
                 fontdict={'fontsize': 22})
    #fig.show()
    filesplit = os.path.splitext(file)
    filename = filesplit[0] + '_map.png'
    fig.subplots_adjust(bottom=0.02, top=0.95, left=0.02, right=1)
    fig.savefig(filename, dpi=600)
    del rast

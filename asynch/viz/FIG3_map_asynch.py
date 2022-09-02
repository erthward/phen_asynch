import matplotlib.pyplot as plt
import matplotlib as mpl
import geopandas as gpd
import rioxarray as rxr
import numpy as np
import palettable
import os, re

countries = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

hotspots = gpd.read_file(('/home/deth/Desktop/CAL/research/projects/'
                          'seasonality/seasonal_asynchrony/analysis/data/hotspots/'
                          'hotspots_2016_1.shp'))

data_dir = '/home/deth/Desktop/CAL/research/projects/seasonality/results/maps/'
files = [f for f in os.listdir(data_dir) if re.search('_global_asynch', f)]

for file in files:

    fig = plt.figure(figsize=(16,8))
    ax = fig.add_subplot(1,1,1)
    rast = rxr.open_rasterio(os.path.join(data_dir, file), masked=True)[2]
    rast = rast.rio.write_crs(4326).rio.reproject(8857)
    cmap = mpl.cm.cubehelix.copy()
    cmap.set_bad('#717171')
    try:
        rast.plot.imshow(ax=ax,
                     zorder=0,
                     cmap=cmap,
                     vmin=np.nanpercentile(rast, 1),
                     vmax=np.nanpercentile(rast, 99),
                     add_colorbar=True,
                    )
    except AttributeError as e:
        print('\n\nAttributeError thrown:\n\t%s\n\n' % e)
    countries.to_crs(8857).plot(ax=ax,
                                color='none',
                                edgecolor='black',
                                linewidth=1,
                                alpha=0.8,
                                zorder=1,
                               )
    #hotspots.to_crs(8857).plot(ax=ax,
    #                          color='none',
    #                          edgecolor='yellow',
    #                          linewidth=1,
    #                          alpha=0.75,
    #                          zorder=2,
    #                         )
    ax.set_xlim(rast.rio.bounds()[0::2])
    ax.set_ylim(rast.rio.bounds()[1::2])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_title('asynch: ' + os.path.splitext(file)[0].split('_')[0],
                 fontdict={'fontsize': 22})
    #fig.show()
    filesplit = os.path.splitext(file)
    filename = filesplit[0] + '_map.png'
    fig.subplots_adjust(bottom=0.02, top=0.95, left=0.02, right=1)

    fig.savefig(filename, dpi=600)
    del rast

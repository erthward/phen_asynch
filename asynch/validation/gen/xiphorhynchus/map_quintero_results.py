import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, MultiPolygon, LineString
import rioxarray as rxr
import alphashape
import matplotlib.pyplot as plt
import sys

sys.path.insert(1, ('/home/deth/Desktop/CAL/research/projects/seasonality/'
                    'seasonal_asynchrony/etc/'))
import phen_helper_fns as phf

# table with geog locs
t2 = pd.read_excel('./Table S2 (Monthly Cloud-Cover).xls')

# table with Mantel results
t3 = pd.read_excel('./Table S3 (Mantel results).xls')

# get all species' number of signif results
ct_signif = t3.groupby('Species').count().loc[:, 'signif'].reset_index()

# join on locs
df = pd.merge(t2.iloc[:, :3], ct_signif, on='Species', how='left')

# get centroids
cent_df = df.loc[:, ['Species', 'Longitude', 'Latitude',
                     'signif']].groupby('Species').mean().reset_index()

# fit hulls
hulls = []
for sp in cent_df['Species'].values:
    try:
        coords = df[df['Species']==sp].loc[:, ['Longitude', 'Latitude']].values
        # calculate alpha of observation coordinates
        hull =  alphashape.alphashape(np.array(coords), alpha=0.25)
        if not isinstance(hull, MultiPolygon):
            if isinstance(hull, LineString):
                hull = hull.buffer(0.05)
            hull = MultiPolygon([hull])
        assert isinstance(hull, MultiPolygon)
    except Exception:
        hull = Point(*cent_df[cent_df['Species'] == sp].loc[:, ['Longitude',
                                                             'Latitude']].values).buffer(0.05)
    hulls.append(hull)

# coerce to geodataframe
cent_df['geom'] = [Point(r['Longitude'],r['Latitude']) for i, r in df.iterrows()]
cent_df['geometry'] = hulls
gdf = gpd.GeoDataFrame(cent_df, geometry='geometry', crs=4326)

# map
asynch = rxr.open_rasterio(phf.ASYNCH_FILES[100])[0]
fig, ax = plt.subplots(1,1)
gdf.plot('signif',
         legend=True,
         cmap='plasma',
         ax=ax,
         edgecolor='black',
         zorder=1,
         alpha=0.5,
        )
xlim = ax.get_xlim()
ylim = ax.get_ylim()
asynch.plot.imshow(ax=ax, cmap='magma', vmin=0, vmax=np.nanpercentile(asynch, 95))
for i, row in gdf.iterrows():
    if pd.notnull(row['signif']):
        ax.text(row['Longitude'],
                row['Latitude'],
                row['Species'],
                size=10,
               )
ax.set_xlim(xlim)
ax.set_ylim(ylim)

fig.savefig('quintero_results_map.png', dpi=700)



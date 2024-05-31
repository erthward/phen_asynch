import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
import h3
import os, re

# NOTE: based on code taken from:
    # https://towardsdatascience.com/uber-h3-for-data-analysis-with-python-1e54acdcc908

# load climate-distance analysis results
res_df = gpd.read_file('./clim_dep_all_MMRR_results_100kmrad.shp')

# load country boundaries
countries = gpd.read_file('../../../data/bounds/NewWorldFile_2020.shp')

# make dataframe to hold h3-converted data
h3_df = pd.DataFrame([], columns=['row_id', 'h3_id',
                                  'h3_geo_boundary', 'h3_centroid'])

# loop over results rows and convert to H3 hexes
for i, row in res_df.iterrows():
    p = row['geometry']
    if isinstance(p, shapely.geometry.MultiPolygon):
        ps = list(p)
    else:
        ps = [p]
    for poly in ps:
        poly_json = gpd.GeoSeries([poly]).__geo_interface__
        poly_json = poly_json['features'][0]['geometry']
        h3_hexes = h3.polyfill_geojson(poly_json, 3)
        for h3_hex in h3_hexes:
            h3_geo_boundary = shapely.geometry.Polygon(
                h3.h3_to_geo_boundary(h3_hex,geo_json=True))
            h3_centroid = h3.h3_to_geo(h3_hex)
            h3_df.loc[len(h3_df)] = [i, h3_hex, h3_geo_boundary, h3_centroid]

# coerce to GeoDataFrame
geoms = [shapely.geometry.Polygon(row['h3_geo_boundary']) for i,
                                            row in h3_df.iterrows()]
h3_df['geometry'] = geoms
h3_gdf = gpd.GeoDataFrame(h3_df, geometry='geometry', crs=4326)

# deduplicate hexes
h3_gdf = h3_gdf.drop_duplicates(subset='geometry')

# summarize results within hexes
mean_clim_dist_vals = []
for i, row in h3_gdf.iterrows():
    mean_clim_dist_val = float(res_df[res_df.intersects(row['geometry'])].mean()['clim_dist'])
    mean_clim_dist_vals.append(mean_clim_dist_val)
assert len(mean_clim_dist_vals) == len(h3_gdf)
h3_gdf['clim_dist_mean'] = mean_clim_dist_vals

# transform to equal-area projection and plot
fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(111)
countries.to_crs(8857).plot(color='none',
                            edgecolor='black',
                            linewidth=1,
                            zorder=0,
                            ax=ax,
                            )
h3_gdf.to_crs(8857).plot('clim_dist_mean',
                         cmap='magma',
                         alpha=0.9,
                         ax=ax,
                         zorder=1)
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_xticks(())
ax.set_yticks(())
ax.set_title('')
fig.subplots_adjust(top=1, bottom=0, left=0, right=1)

fig.savefig('clim_dep_hex_summary.png', dpi=700)

import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import numpy as np
import os

###############################################################################
# NOTE: MANUAL STEPS TAKEN:
    # 1. capture screenshot of Fig. 1a in Zhisheng et al. 2015, "Global Monsoon
    #    Dynamics and Climate Change"
    # 2. go to https://apps.automeris.io/wpd/, click File -> Load Images and upload screenshot
    # 3. select "Map With Scale Bar", click "Align Axes", then drop two points
    #    along the 0-degree meridian (points: (0,0) and (0, 30)

    # 4. click "Complete" and enter '10' and 'degrees' in the prompt boxes

    # 5. click on the color box next to "Foreground color" and select the blue
    #    corresponding to the DJF ITCZ curve

    # 6. set color distance to 20, deltaX and deltaY to 1 pixel, then click
    #    "Run" to add automatically digitized points along the DJF ITCZ curve

    # 7. click "Delete Point (D)", hover over the incorrectly digitized points
    #    added to the map legend, and repeatedly click until the last point
    #    there is deleted

    # 8. click "View Data" -> "Download .CSV" and download the CSV file as
    #    'ITCZ_DJF.csv'

    # 9. repeat steps 5-8 for the JJA ITCZ curve and download its file as
    #    'ITCZ_JJA.csv'

    # 10. run this script to produce the shapefile

###############################################################################

data_dir = ('/home/deth/Desktop/CAL/research/projects/seasonality/'
            'seasonal_asynchrony/analysis')

# read in raw, digitized data for boreal-summer ITCZ
itczJJA = pd.read_csv(os.path.join(data_dir, 'ITCZ_JJA.csv'), header=None)
itczJJA.columns = ['lon', 'lat']
# drop the few points that were digitized to the legend line
# (rather than the actual mapped line)
itczJJA = itczJJA[itczJJA['lat']<140]
# and average all points that wound up having identical x values
itczJJA = itczJJA.groupby('lon').mean().reset_index()
# add minute jitter to lons, so that none are identical
#itczJJA['lon'] = itczJJA['lon'] + np.random.normal(0, 0.0000001,
#                                                  size=len(itczJJA))
#itczJJA = itczJJA.sort_values('lon')

# same for boreal-winter ITCZ
itczDJF = pd.read_csv(os.path.join(data_dir, 'ITCZ_DJF.csv'), header=None)
itczDJF.columns = ['lon', 'lat']
itczDJF = itczDJF[itczDJF['lat']<140]
itczDJF = itczDJF.groupby('lon').mean().reset_index()
#itczDJF['lon'] = itczDJF['lon'] + np.random.normal(0, 0.0000001,
#                                                   size=len(itczDJF))
#itczDJF = itczDJF.sort_values('lon')

# use the fact that the JJA line first and last points are at 0deg latitude,
# and the fact that points were digitized as increasing downward, so need to be
# 'flipped', to correct all points' latitudes
JJA_DJF_offset = itczJJA['lat'][0] - itczDJF['lat'][0]
itczJJA['lat'] = itczJJA['lat'][0] - itczJJA['lat']
itczDJF['lat'] = itczDJF['lat'][0] - itczDJF['lat'] + JJA_DJF_offset

# rescale lon values to be between -180 and 180
rescaler = MinMaxScaler((-180, 180))
itczJJA['lon'] = rescaler.fit_transform(itczJJA['lon'].values.reshape(-1, 1))
itczDJF['lon'] = rescaler.fit_transform(itczDJF['lon'].values.reshape(-1, 1))

# combine into geodataframe
itczJJA_line = LineString(itczJJA.values)
itczDJF_line = LineString(itczDJF.values)
df = {'time': ['JJA', 'DJF'], 'geometry': [itczJJA_line, itczDJF_line]}
gdf = gpd.GeoDataFrame(df, crs='EPSG:4326')

# write out
gdf.to_file('ITCZ_li_zeng_2005_digitized.shp')

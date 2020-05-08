#!/usr/bin/python

# from:
    # https://www.tylanderson.com/post/uploading-custom-raster-image-collections-to-google-earth-engine/

from osgeo import gdal
import pandas as pd
import os

in_dir = (r"/run/media/drew/SLAB/seasonality/SIF/OCO-2/"
          "gridded/Global_High_Res_SIF_OCO2_1696/data")
out_csv = r"./SIF_OCO2_ANN-gridded_metadata.csv"

# List of metadata keys to extract
meta_keys = ['cloud_coverage',
             'SENSOR',
             'spatial_coverage']

# Add other required keys and make pandas data frame
# We will store the metadata in the df
csv_keys = ['id_no', 'system:time_start']
csv_keys.extend(meta_keys)
df = pd.DataFrame(columns=csv_keys)

for root, dirs, files in os.walk(in_dir):
    for file in files:
        if file.endswith(".tif"):
            path = os.path.join(root, file) # make path
            raster = gdal.Open(path, gdal.GA_ReadOnly) # Open Raster in GDAL
            meta = dict.fromkeys(csv_keys) # Make a dict to store key, values

            # Get time from filename
            # Example filename (LC80300062013106C2V01_MTLstack.tif)
            year = file[9:13]
            doy = file[13:16]
            hour = '12'
            timestamp = year + doy + hour
            # convert time stamp to gee format
            gee_timestamp = int(pd.to_datetime(timestamp,
                                               format='%Y%j%H').value // 10**6)

            meta['id_no'] = os.path.splitext(file)[0]
            meta['system:time_start'] = gee_timestamp
            for key in meta_keys:
                meta[key] = raster.GetMetadataItem(key)

            df = df.append(meta, ignore_index=True)
print(df)
df.to_csv(out_csv, index=False)

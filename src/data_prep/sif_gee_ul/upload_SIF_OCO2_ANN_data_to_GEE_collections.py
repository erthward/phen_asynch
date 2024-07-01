#!/usr/bin/python

import os
import sys
import glob
import pandas as pd
import rasterio as rio

dataset = sys.argv[1].lower()
assert dataset in ['sif'], 'Valid datasets include: "SIF".'

# Adapted from:
    # https://www.tucson.ars.ag.gov/notebooks/uploading_data_2_gee.html

# Define static variables
gc_bucket_dict = {'sif': "seasonality_data/OCO2_SIF_ANN",
                 }
gc_bucket = gc_bucket_dict[dataset]
print(gc_bucket)

img_coll_dict = {'sif': "users/drewhart/seasonality_data/OCO2_SIF_ANN",
                }
img_coll = img_coll_dict[dataset]

# Create the ImageCollection asset in GEE
#cmd = 'earthengine create collection %s' % img_coll
#print("\n Creating ImageCollection using the following command:\n\t%s\n" % cmd)
#os.system(cmd)

# Get the provider and URL strings
str_prov_dict = {'sif': ("(string)provider=Oak Ridge National "
                         "Laboratory (ORNL) Distributed Active "
                         "Archive Center (DAAC)"),
                }
str_prov = str_prov_dict[dataset]

url_dict = {'sif': ('(string)URL=https://daac.ornl.gov/VEGETATION/guides/'
                    'Global_High_Res_SIF_OCO2.html'),
           }
url = url_dict[dataset]

# set the data directory
data_dir_dict = {'sif': ("/run/media/drew/SLAB/seasonality/SIF/OCO-2/"
                         "gridded/Global_High_Res_SIF_OCO2_1696/data"),
                }
data_dir = data_dir_dict[dataset]

# read in the DataFrame of file start and end dates
df_file_dict = {'sif': 'SIF_OCO2_ANN_upload_metadata.csv',
               }
df_file = df_file_dict[dataset]
df = pd.read_csv(os.path.join(data_dir, df_file))

# set the CRS
crs_dict = {'sif': 'EPSG:4326',
           }
crs = crs_dict[dataset]

# Get file names to extract date and call ingestion command for each file to be added into an asset as image collection
# Example of filenames used here are from a monthly timeseries: rainfall_US_20140101.tif
for f in glob.glob(os.path.join(data_dir, '*.tif')):
    # get the nodata val
    ds = rio.open(f)
    nodata_val = ds.nodata

    # get file basename and its row of data
    basename = os.path.splitext(os.path.split(f)[-1])[0]
    row = df[df['id_no'] == basename]
    print(f'NOW MOVING: {basename}...\n')

    # get the start and end timestamps
    start = str(row['system:time_start'].values[0])
    end = str(row['system:time_end'].values[0])

    # create asset name
    asset=os.path.join(img_coll, basename)

    # build the earthengine command
    cmd = (f'earthengine upload image --asset_id="{asset}" '
           f'--pyramiding_policy=mean --time_start="{start}" --time_end="{end}"'
           f" --property='{str_prov}' --property='{url}' "
           f'--nodata_value={nodata_val} '
           f'--crs={crs} gs://{os.path.join(gc_bucket, os.path.split(f)[-1])}')
    test_cmd = 'echo ' + cmd
    os.system(test_cmd)
    # Call the ingestion command to create the Asset and add to the
    # corresponding ImageCollection
    os.system(cmd)
    print('Command:\n--------')
    print(f"{'='*80}\n")

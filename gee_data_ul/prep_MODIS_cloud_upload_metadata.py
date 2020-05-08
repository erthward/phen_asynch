#!/usr/bin/python

import pandas as pd
import glob
import os

# define directory, files, output file, output DataFrame
data_dir = "/run/media/drew/SLAB/seasonality/other/cloud"
files = glob.glob(os.path.join(data_dir, '*.tif'))
print(files)

output_file = './MODIS_cloud_upload_metadata.csv'
output_cols = ['id_no', 'system:time_start', 'system:time_end']
output_df = pd.DataFrame(columns=output_cols)

gee_start_timestamp = '2000-01-01T00:00:00'
gee_end_timestamp = '2015-01-01T00:00:00'
#gee_start_date = '2000-01-01 00:00:00'
#gee_end_date = '2015-01-01 00:00:00'
# NOTE: geeup wound up not working on my machine, so no need for the
# milliseconds-since-UNIX timestamp; plain dates will do just fine
#gee_start_timestamp = pd.to_datetime(gee_start_date,
#                                     format='%Y-%m-%d %H:%M:%S').value // 10**6
#gee_end_timestamp = pd.to_datetime(gee_end_date,
#                                   format='%Y-%m-%d %H:%M:%S').value // 10**6
# NOTE: For some reason, these dates are coming out 8 hours too early;
# not sure why, but for now I'll just add 8 hours' worth of seconds onto them
#gee_start_timestamp += 8*60*60*1000
#gee_end_timestamp += 8*60*60*1000

for f in files:
    # create a row of data for this file
    row = dict.fromkeys(output_cols)
    row['id_no'] = os.path.splitext(os.path.split(f)[-1])[0]
    row['system:time_start'] = gee_start_timestamp
    row['system:time_end'] = gee_end_timestamp
    # append the row to the output CSV
    output_df = output_df.append(row, ignore_index=True)

print(output_df)

# write the df to disk
output_df.to_csv(output_file, index=False)

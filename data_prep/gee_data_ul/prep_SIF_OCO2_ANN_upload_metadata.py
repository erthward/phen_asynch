#!/usr/bin/python

import pandas as pd
import datetime
import glob
import os
import re

# define directory, files, output file, output DataFrame
data_dir = ("/run/media/drew/SLAB/seasonality/SIF/OCO-2/"
            "gridded/Global_High_Res_SIF_OCO2_1696/data")
files = glob.glob(os.path.join(data_dir, '*.nc'))
#files = glob.glob(os.path.join(data_dir, '*.tif'))
print(files)

output_file = './SIF_OCO2_ANN_upload_metadata.csv'
output_cols = ['id_no', 'system:time_start', 'system:time_end']
output_df = pd.DataFrame(columns=output_cols)


# function to get a file's date and return it as a string with
# the correct date formatting
def get_file_gee_timestamp(f, start=True):
    # regex for the date in the filename
    date_patt = '(?<=ann_)\d{6}[ab]'
    date_str = re.search(date_patt, f).group()
    # format as a string, making the day 1 or 16, depending on 'a' or 'b' file
    date = '%s-%s-%s' % (date_str[:4],
                         date_str[4:6],
                         str((1 + (15 * (date_str[-1] == 'b')))).zfill(2))
    # convert to the end-time, if necessary
    if not start:
        if date_str[-1] == 'a':
            before_len = len(date)
            date = date[:8] + '16'
            after_len = len(date)
            assert before_len == after_len
        else:
            if int(date[5:7]) < 12:
                before_len = len(date)
                date = date[:5] + str(int(date[5:7]) + 1).zfill(2) + '-01'
                after_len = len(date)
                assert before_len == after_len
            else:
                before_len = len(date)
                date = str(int(date[:4]) + 1) + '-01-01'
                after_len = len(date)
                assert before_len == after_len
    # turn into date object in GEE format
    date = date + 'T00:00:00'

    # NOTE: geeup wound up not working on my machine, so no need for the
    # milliseconds-since-UNIX timestamp; just plain dates will do
    #gee_timestamp = pd.to_datetime(date,
    #                               format='%Y-%m-%d %H:%M:%S').value // 10**6
    # NOTE: for some reason all my timestamps were coming out as 17:00:00 hours;
    # don't know why, but for now I'm just subtracting 17 hours' worth of
    # milliseconds
    #gee_timestamp = gee_timestamp - (17*60*60*1000)

    #return gee_timestamp
    return date


for f in files:
    # create a row of data for this file
    row = dict.fromkeys(output_cols)
    row['id_no'] = os.path.splitext(os.path.split(f)[-1])[0]
    row['system:time_start'] = get_file_gee_timestamp(f)
    row['system:time_end'] = get_file_gee_timestamp(f, start=False)
    # append the row to the output CSV
    output_df = output_df.append(row, ignore_index=True)

print(output_df)

# write the df to disk
output_df.to_csv(output_file, index=False)

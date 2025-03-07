import os
import re
import json
import urllib.request
import numpy as np
import pandas as pd

# set data dir
data_dir = '/media/deth/SLAB/diss/3-phn/phenocam/'

# download data?
download = False

# read all the ROIs names
rois = pd.read_csv('./ROIs.csv')
roi_names = np.unique(rois['roi_name'].values)

# dict to store info from meta file
meta_info = {'sitename': [],
             'long_name': [],
             'lat': [],
             'lon': [],
             'elevation': [],
             'date_start': [],
             'date_end': [],
             'MAP_worldclim': [],
             'MAT_worldclim': [],
             'primary_veg_type': [],
             'secondary_veg_type': [],
             'landcover_igbp': [],
             'infrared_enabled': [],
             'site_acknowledgements': [],
            }

# try to grab an NDVI 3-day file for each one,
# skipping errors
for roi_name in roi_names:
    try:
        url = f"https://phenocam.nau.edu/data/archive/{roi_name.split('_')[0]}/ROI/"
        filename = f'{roi_name}_ndvi_3day.csv'
        fileurl = os.path.join(url, filename)
        filepath = os.path.join(data_dir, filename)
        # NOTE: -c to clobber, -P to specify download destination dir
        wget_cmd = f"wget -c -P {data_dir} {fileurl}"
        if download:
            os.system(wget_cmd)
        # also populate the metadata table that we'll need to add to supps
        if roi_name.split('_')[0] not in meta_info['sitename']:
            meta_url = f"https://phenocam.nau.edu/data/archive/{roi_name.split('_')[0]}/"
            meta_filename = f"{roi_name.split('_')[0]}_meta.json"
            meta_fileurl = os.path.join(meta_url, meta_filename)
            with urllib.request.urlopen(meta_fileurl) as meta:
                meta_dict = json.load(meta)
                for k in meta_info:
                    try:
                        val = meta_dict[k]
                    except Exception:
                        val = meta_dict['phenocam_site'][k]
                    meta_info[k].append(val)
    except Exception as e:
        print(e)
        for k in meta_info:
            if k == 'sitename':
                meta_info[k].append(roi_name.split('_')[0])
            else:
                meta_info[k].append(np.nan)

# save metadata table to file
meta_df = pd.DataFrame(meta_info)
meta_df.to_csv('site_metadata.csv', index=False)

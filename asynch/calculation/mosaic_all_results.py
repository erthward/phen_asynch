# mosaic_all_results.py

import os

overwrite = True

# map layers' arg strings to their filename strings
lyrs_dict = {'c': 'coeffs',
             'r2': 'harm_reg_R2',
             'a': 'asynch',
            }
# map neigh_rads' arg strings to their filename strings
neigh_rads_dict = {str(rad): '_%ikm' % rad for rad in [50, 100, 150]}
neigh_rads_dict[''] = ''

# template command for mosaicking script
template_cmd = "python %s %s %s %s %s"

# get path for savio...
if os.getcwd().split('/')[1] == 'global':
    script_path = ("/global/home/users/drewhart/seasonality/seasonal_asynchrony/"
                   "asynch/calculation/mosaic_results.py")
    data_path = "/global/scratch/users/drewhart/seasonality/GEE_outputs/"
# ... or on laptop
else:
    script_path = ("/home/deth/Desktop/CAL/research/projects/seasonality/"
                   "seasonal_asynchrony/asynch/calculation/mosaic_results.py")
    data_path = "/media/deth/SLAB/diss/3-phn/GEE_outputs/"


# loop over vars and lyrs
for var in ['NIRv', 'NIRv_STRICT', 'SIF', 'SIF_STRICT',
            'tmmn', 'pr', 'def', 'cloud']:

    print(f"MOSAICKING ALL DATA FOR {var}...\n\n")

    # get this var's specific data path
    data_subpath = os.path.join(data_path, var)
    if len(os.listdir(data_subpath)) > 0:

        for lyr in ['c', 'r2', 'a']:

            # determine neighborhood-radius strings to be used, depending on lyr
            if lyr == 'a':
                neigh_rads = ['50', '100', '150']
            else:
                neigh_rads = ['']

            for neigh_rad in neigh_rads:
                try:
                    # set output filepath
                    filename = f"{var}_{lyrs_dict[lyr]}{neigh_rads_dict[neigh_rad]}.tif"
                    filepath = os.path.join(data_subpath, filename)

                    if (overwrite or
                        (not overwrite and
                         not os.path.isfile(filepath))):

                        # format the command string
                        cmd_str = template_cmd % (script_path,
                                                  data_subpath,
                                                  filepath,
                                                  lyr,
                                                  neigh_rad,
                                                 )
                    else:
                        print((f"\t{filepath} ALREADY EXISTS AND SCRIPT RUN WITH "
                               "overwrite=False. SKIPPING.\n\n"))
                    # call command
                    print("\n\nNOW CALLING:\n\t%s\n\n" % cmd_str)
                    os.system(cmd_str)

                except Exception as e:
                    print("\n\nCOULD NOT MOSAIC FOR %s, %s %s. MOVING ON...\n\n" % (
                        var,
                        lyr,
                        (('(neigh_rad=%s)' % neigh_rad) * (neigh_rad != ''))))

            print('\n\n~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~.~\n')

    else:
        print(f"\tNO DATA FOUND FOR {var}. SKIPPING.\n\n")

print("\n\nALL RESULTS MOSAICKED SUCCESSFULLY!\n\n")

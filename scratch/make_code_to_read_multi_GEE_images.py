import numpy as np
import os
import sys

# frequency (i.e. every 'n' days) for which to produce code to read
# set of annual images
freq = int(sys.argv[1])

# main command
main_cmd = 'var ic = ee.ImageCollection([\n\t%s]);'

# generate all Image-reading commands
days = np.arange(0, 365+freq, freq)
if days[-1] > 364:
    if days[-2] == 364:
        days = days[:-1]
    else:
        days[-1] = 364
image_read_cmds = ['ee.Image("users/ldimaggio/SIF%i")' % day for day in days]

# concatenate all cmds and insert into main command
full_cmd = main_cmd % ',\n'.join(image_read_cmds)

# write to file
with open('GEE_cmd_every_%i_days.txt' % freq, 'w') as f:
    f.write(full_cmd)

print('Booyah.')

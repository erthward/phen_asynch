#!/usr/bin/python
#download_data.py

'''Download all the OCO-2 SIF data I need.'''

import os

#grab pwd
with open('/home/ihavehands/Desktop/NASA_pw.txt', 'r') as f:
    pwd = f.read().rstrip('\n')

#set the NASA data-archive URL
url = 'https://oco2.gesdisc.eosdis.nasa.gov/data/s4pa/OCO2_DATA/OCO2_L2_Lite_SIF.8r/'

#create the wget command, into which the auth_opts, pw, and formatted archive url will be added
cmd_str = 'wget --content-disposition %s --user=drew.hart --password=%s -r -nd -c -nH -np -A nc4,xml %s'

#set authorization options
auth_opts = '--load-cookies ~/Desktop/urs_cookies --save-cookies ~/Desktop/urs_cookies --keep-session-cookies'

#create the full command
cmd = cmd_str % (auth_opts, pwd, url)

#print the command
print('PRESS <Enter> TO RUN THE FOLLOWING COMMAND:\n\t%s\n\n' % cmd)
input()

#run the command
os.system(cmd)

print('\n\nFINISHED DOWNLOADING ALL DATA.\n\n\n')


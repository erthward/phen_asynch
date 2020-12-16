import numpy as np
import asynch_fns as af
import tensorflow as tf
import os

kernel_size = 60

os.chdir('/home/drew/Desktop/TEST_upload_asynch')


def drop_kernel_overlap(filename, kernel_size):
    print('\nnow dropping kernel overlap for file %s\n' % filename)
    kernel_hw = int(kernel_size/2)
    ip, op = af.get_inpatches_outpatches(filename, af.OUTBANDS, af.DIMS)
    # drop the overlap margin
    patches = [p[:,
                 kernel_hw:-kernel_hw,
                 kernel_hw:-kernel_hw,] for p in ip]
    out_filename = os.path.splitext(filename)[0] + '_DROP_OVERLAP.tfrecord'
    af.write_tfrecord_file(patches, out_filename, af.OUTBANDS)


for filename in os.listdir('.'):
    if (not 'mixer' in filename and
        filename.endswith('tfrecord') and
        not 'DROP_OVERLAP' in filename):
        drop_kernel_overlap(filename, kernel_size)

print('DONE')

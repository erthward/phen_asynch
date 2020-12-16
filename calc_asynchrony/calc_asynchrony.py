#!/bin/python
# calc_asynchrony.py

# TODO:
    # 0. figure out how to scale up permutation tests on GEE
    # 0.5 chat with Ian about asynch, R2s, ns, and scaling limitations
    # 1. figure out what to do about negative R2 values
    # 2. figure out how to write out the TFRecord files
    # 3. begin planning how to parallelize on XSEDE/Savio
    # 4  figure out how to export the coeff patches globally, to Google Drive,
    #    then copy over to XSEDE/Savio from there

#--------
# imports
#--------

# import our module full of functions
import asynch_fns as af

# import other packages
import multiprocessing as mp
import sys

# limit number of files to process in this job, if necessary
if len(sys.argv) > 1:
    try:
        start_filenum_idx = int(sys.argv[1])
        stop_filenum_idx = int(sys.argv[2])
        stop_filenum_idx = max(stop_filenum_idx, len(af.FILENAMES)-1)
    except Exception as e:
        raise ValueError(('If providing starting and ending file-number indices'
                          ' then they must be provided as integers as the two '
                          'arguments immediately following the script\'s '
                          'filename.'))
else:
    start_filenum_idx = 0
    stop_filenum_idx = len(af.FILENAMES)-1

if __name__ == '__main__':

    # get the files_dict to process for this job
    # (using starting and ending files' indices)
    files_dict = {fn: af.FILES_DICT[fn] for fn in af.FILENAMES[
                                            start_filenum_idx:stop_filenum_idx]}
    print('\n\nFiles to be processed:\n\n\t - %s\n\n' % ('\n\n\t - '.join(
                                                    [fn for fn in files_dict])))

    print('-'*80, '\n\n')
    # how many CPUs?
    ncpu = mp.cpu_count()
    print('%i CPUs AVAILABLE IN TOTAL\n\n' % ncpu)
    ncpu = min(ncpu, len(files_dict))
    print('THIS COMPUTATION WILL BE COMPLETED USING %i CPUS\n\n' % ncpu)
    print('-'*80, '\n\n')

    # set the start method to 'spawn' instead of 'fork, to avoid deadlock
    # (see: https://pythonspeed.com/articles/python-multiprocessing/)
    mp.set_start_method('spawn')

    # make the Pool
    pool = mp.Pool(ncpu)

        # map the input files dict into the CPUs in our pool
    print("BEGIN ASSIGNING JOBS...\n")
    pool.map_async(af.main_fn, files_dict.items())

    # close and join the Pool
    pool.close()
    pool.join()

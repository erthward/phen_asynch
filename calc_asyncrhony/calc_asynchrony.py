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

if __name__ == '__main__':

    # how many CPUs?
    ncpu = mp.cpu_count()
    ncpu = min(ncpu, len(af.FILES_DICT))
    print('\n\nTHIS COMPUTATION WILL BE COMPLETED USING %i CPUS\n\n' % ncpu)

    # set the start method to 'spawn' instead of 'fork, to avoid deadlock
    # (see: https://pythonspeed.com/articles/python-multiprocessing/)
    mp.set_start_method('spawn')

    # make the Pool
    #pool = pomp.ProcessingPool(ncpu)
    pool = mp.Pool(ncpu)
    
    # map the input files dict into the CPUs in our pool
    print("BEGIN ASSIGNING JOBS...\n")
    #pool.amap(main_fn, files_dict.items())
    #jobs = []
    #for item in files_dict.items():
    #    job = pool.apply_async(main_fn, args=[item])
    #    jobs.append(job)
    pool.map_async(af.main_fn, af.FILES_DICT.items())
    
    # close and join the Pool
    pool.close()
    pool.join()
    

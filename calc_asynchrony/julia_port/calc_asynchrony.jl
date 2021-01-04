"""
NOTE: should be run after calling `julia -p <NCPU-1>`

Loads the asynchrony functions and constants into all available workers'
workspaces, then maps the main_fn across those workers.

This will calculate asynchrony and write output files for all available input TFRecord files.
"""

# load the Distributed package
using Distributed

# include the asynch functions
@everywhere include("./asynch_fns.jl")

#--------------------------------------------------
# use parallelization to call main_fn for all files
#--------------------------------------------------
println("CPU INFO:\n--------")
println(Sys.cpu_info())
println("\n\n")
np = nprocs()
nw = nworkers()
println("USING $np PROCESSES, WITH $nw WORKERS")

pmap(main_fn, zip(keys(FILES_DICT), values(FILES_DICT)))

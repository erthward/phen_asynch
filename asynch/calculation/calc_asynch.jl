"""
NOTE: should be run after calling `julia -p <NCPU-1>`

Loads the asynchrony functions and constants into all available workers'
workspaces, then maps the main_fn across those workers.

This will calculate asynchrony and write output files for all available input TFRecord files.
"""

# load the Distributed package
using Distributed
using ArgParse

# include the asynch functions
if splitpath(pwd())[3] == "home" || splitpath(pwd())[3] == "scratch"
    @everywhere include("/global/home/users/drewhart/seasonality/seasonal_asynchrony/" *
		        "asynch/calculation/asynch_fns.jl")
else
    @everywhere include("/home/deth/Desktop/CAL/research/projects/seasonality/" *
                        "seasonal_asynchrony/asynch/calculation/asynch_fns.jl")
end


# parse the command-line args
function parse_cmdline_args()
    s = ArgParseSettings()

    @add_arg_table! s begin
        # get the variable for which to calculate asynchrony
        "--var"
            help = """which variable's files to process

                      valid values: NIRv, NIRv_STRICT, SIF, SIF_STRICT, tmmn, pr, def, cloud
                   """
            arg_type = String
            default = "NIRv"
        # get the neighborhood radius within which to calculate asynchrony
        "--neigh_rad"
            help = """neighborhood radius (in km) to use

                      valid values: 50, 100, 150
                   """
            arg_type = String
            default = "100"

    end

    return parse_args(s)
end


#--------------------------------------------------
# use parallelization to call main_fn for all files
#--------------------------------------------------
function main()

    # parse the args
    parsed_args = parse_cmdline_args()

    # process arg values
    var = parsed_args["var"]
    @assert var in ["NIRv" "NIRv_STRICT" "SIF" "SIF_STRICT" "def" "pr" "tmmn" "tmmx" "cloud"]
    neigh_rad = parsed_args["neigh_rad"]
    @assert neigh_rad in ["50", "100", "150"]
    # convert to meters
    neigh_rad = 1000 * parse(Int64, neigh_rad)
    println("\nCALCULATING FOR $var WITHIN $neigh_rad-METER NEIGHBORHOOD...\n")
   
    # use var to get data dir
    abs_data_dir = BASE_DATA_DIR * var
    # NOTE: the TFRecord package throws a globbing error
    # ("ERROR: Glob pattern cannot be empty or start with a / character")
    # when I feed TFRecordReader an absolute path
    # (SEEMS LIKE A BUG, NO?!),
    # so get the relative path to the data_dir instead
    data_dir = relpath(abs_data_dir)

    # get files' row, col, and patch info as a Dict containing input filenames
    # as keys and the files' patch, column, and row numbers as values.
    # Needed in order to parallelize the computation across files
    # while still calculating lats and long correctly for each file.
    FILES_DICT = get_row_col_patch_ns_allfiles(data_dir, PATT_AFT_FILENUM, neigh_rad)
    files_dict = deepcopy(FILES_DICT)

    @info "\n$(length(files_dict)) FILES TO BE PROCESSED:\n----------------------\n"
    for (k,v) in sort(files_dict)
        # NOTE: make sure we did not grab any previous output files
        @assert ! occursin("-OUT-", k)
        println("\t$k")
    end

    # display CPU and parallelization info
    @info "CPU SUMMARY:\n--------" 
    println(Sys.cpu_info())
    println("\n\n")
    np = nprocs()
    nw = nworkers()
    @info "USING $np PROCESS$(np > 1 ? "ES" : ""), WITH $nw WORKER$(nw > 1 ? "S" : "")"
       
    pmap(main_fn, zip(repeat([data_dir], length(files_dict)),
                      repeat([neigh_rad], length(files_dict)),
                      keys(files_dict),
                      values(files_dict)))
end

main()

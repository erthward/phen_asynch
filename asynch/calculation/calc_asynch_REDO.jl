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
        # TODO: ADD RAD AND VAR CMDLINE ARGS HERE
        "--which"
            help = """which files to process

                      'all' will process everything

                      'unprocessed' will process files without output files,
                                    or files with output files of size 0b

                      TODO: add handling for number ranges
                   """
            arg_type = String
            default = "all"
    end

    return parse_args(s)
end


#--------------------------------------------------
# use parallelization to call main_fn for all files
#--------------------------------------------------
function main()

    # parse the args
    parsed_args = parse_cmdline_args()

    # use args to determine behavior
    if parsed_args["which"] == "all"
        files_dict = deepcopy(FILES_DICT)
    elseif parsed_args["which"] == "unprocessed"
        files_dict = Dict()
        for (k, v) in FILES_DICT
            outfile_name = replace(k, r"(?<=\d{5})." => "_OUT.")
            if !isfile(outfile_name) || filesize(outfile_name) == 0
                files_dict[k] = v
            end
        end
    end

    #TODO @assert VAR in ["NIRv" "NIRv_STRICT" "SIF" "SIF_STRICT" "def" "pr" "tmmn" "cloud"]
  
    # use var to get data dir
    ABS_DATA_DIR = DATA_DIR * VAR
    # NOTE: the TFRecord package throws a globbing error
    # ("ERROR: Glob pattern cannot be empty or start with a / character")
    # when I feed TFRecordReader an absolute path
    # (SEEMS LIKE A BUG, NO?!),
    # so get the relative path to the data_dir instead
    DATA_DIR = relpath(ABS_DATA_DIR)


    # TODO: MAKE FILES_DICT HERE USING RAD AND VAR ARGS
    # get files' row, col, and patch info as a Dict containing input filenames
    # as keys and the files' patch, column, and row numbers as values.
    # Needed in order to parallelize the computation across files
    # while still calculating lats and long correctly for each file.
    FILES_DICT = get_row_col_patch_ns_allfiles(DATA_DIR, PATT_AFT_FILENUM)





    @info "\n$(length(files_dict)) FILES TO BE PROCESSED:\n----------------------\n"
    for (k,v) in sort(files_dict)
        println("\t$k")
    end

    # display CPU and parallelization info
    @info "CPU SUMMARY:\n--------" 
    println(Sys.cpu_info())
    println("\n\n")
    np = nprocs()
    nw = nworkers()
    @info "USING $np PROCESS$(np > 1 ? "ES" : ""), WITH $nw WORKER$(nw > 1 ? "S" : "")"
       
    #TODO: use repeat([var], length(files_dict)) and same for rad to add them into mainfn as args
    pmap(main_fn, zip(keys(files_dict), values(files_dict)))
end

# TODO: loop main over all vars and rads
# TODO: NEED TO REQUEST EACH LOOP TO WAIT UNTIL COMPLETE?
# TODO: START WITH JUST 2 RADS AND 2 VARS, CLEAR ALL PREV OUTPUT, THEN MOSAIC AND CONFIRM RESULTS BEFORE SCALING UP TO ALL
main()

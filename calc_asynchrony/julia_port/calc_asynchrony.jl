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
@everywhere include("./asynch_fns.jl")

# parse the command-line args
function parse_args()
    s = ArgParseSettings()

    @add_arg_table s begin
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
    parsed_args = parse_args()

    # use args to determine behavior
    if parsed_args["which"] == "all"
        files_dict = deepcopy(FILES_DICT)
    else if parsed_args["which"] == "unprocessed"
        files_dict = Dict()
        for (k, v) in FILES_DICT
            outfile_name = replace(k, r"-(?=\d{5})" => "-OUT-")
            if !isfile(outfile_name) or filesize(outfile_name) == 0
                files_dict[k] = v
            end
        end
    end
    println("\nFILES TO BE PROCESSED:\n----------------------\n")
    for (k,v) in files_dict
        println("\t$k")
    end

    # display CPU and parallelization info
    println("CPU INFO:\n--------")
    println(Sys.cpu_info())
    println("\n\n")
    np = nprocs()
    nw = nworkers()
    println("USING $np PROCESSES, WITH $nw WORKERS")

    pmap(main_fn, zip(keys(files_dict), values(files_dict)))
end

main()

using OrderedCollections
using NearestNeighbors
using StaticArrays
using Distributed
using Statistics
using Distances
using GeoArrays
using TFRecord
using Printf
using Images
using Colors
using Plots
using Dates
using JSON
using Glob
using GLM


#-----------
# set params
#-----------

# stdout options
"""
Whether to use verbose output
"""
const VERBOSE = true
"""
Whether to time the asynchrony calculation
(Time will be displayed per pixel)
"""
const TIMEIT = true


# choose the plots backend
# (sticking with the python plotting window that I know well, for now)
#pyplot()
# NOTE: pyplot breaks when trying to plot the Grayscale image,
# so for now using GR instead
gr()


# the TFRecord package throws a globbing error
# ("ERROR: Glob pattern cannot be empty or start with a / character")
# when I feed TFRecordReader an absolute path
# (SEEMS LIKE A BUG, NO?!),
# so get the relative path to the data_dir instead

"""
absolute path to the data directory
"""
const ABS_DATA_DIR = "/home/drew/Desktop/tmp_output"

"""
pattern that occurs just before the file number
"""
const PATT_B4_FILENUM = "SIF_OUT-"

"""
directory where the TFRecord data and mixerfile live
"""
const DATA_DIR = relpath(ABS_DATA_DIR)

"""
default missing-data val
"""
const DEFAULT_VAL = -9999.0

"""
names of the bands to save in the output TFRecord files
"""
const BANDS = ["asynch", "asynch_R2", "asynch_euc", "asynch_euc_R2", "asynch_n"]


#-----------------
# define functions
#-----------------
#
"""
Uses the data-directory path provided to get and return
lists of the input and output files' paths.
"""
function get_file_paths(data_dir::String)::Array{String,1}
    # set up IO paths
    filepaths = glob("*tfrecord", relpath(data_dir))
    # order the infilepaths
    sort!(filepaths)
    return filepaths
end


"""
Gets the info out of the mixerfile located in the data_dir.
"""
function read_mixer_file(data_dir::String)::Dict{String,Any}
    # get possible mixer filepaths
    mixerfilepaths = glob("*mixer.json", data_dir)
    # make sure there's only 1
    @assert(length(mixerfilepaths) == 1, "Did not find exactly one valid mixerfile!")
    # grab the correct filepath
    mixerfilepath = mixerfilepaths[1]
    # load it into a Dict
    mixer = JSON.Parser.parsefile(mixerfilepath)
    return mixer
end


"""
Parses the needed information out of a Dict of mixerfile contents,
altering it as necessary.
"""
function get_mixer_info(mixer_content::Dict{String, Any})::Tuple{Tuple{Int64, Int64},
                                                                 String,
                                                                 Float64, Float64,
                                                                 Float64, Float64,
                                                                 Int64, Int64}
    # get patch dimensions, adding the kernel size to correct for
    # the error in the GEE TFRecord output
    dims = calc_patch_dimensions(mixer_content)
    # get the SRS and its affine projection matrix
    srs = mixer_content["projection"]
    crs = srs["crs"]
    affine = srs["affine"]["doubleMatrix"]
    # get the x and y resolutions
    xres, yres = affine[1], affine[5]
    # get the x and y min values (i.e. the center coordinates of
    # the top-left pix)
    # NOTE: need to subtract half-kernel-size pixels' worth of res, to account for
    #       fact that mixer file doesn't include the size of the kernel
    #       used to create neighbor-patch overlap
    # NOTE: need to add half pixel's worth of res
    #       (i.e. subtract it from the amount being subtracted)
    #       to account for the fact that the xmin and ymin are expressed
    #       as upper-left pixel corner, whereas I want to operate
    #       on basis of pixel centers
    xmin = affine[3]
    ymin = affine[6]
    # get the number of patches per row and the total number of rows
    # NOTE: FOR NOW, ASSUMING THAT THE MIXER IS SET UP SUCH THAT THE MOSAIC SHOULD
    # BE FILLED ROW BY ROW, LEFT TO RIGHT
    patches_per_row = mixer_content["patchesPerRow"]
    tot_patches = mixer_content["totalPatches"]
    return dims, crs, xmin, ymin, xres, yres, patches_per_row, tot_patches
end


"""
Calculates and returns the patch dimensions,
using the input mixerfile info and the kernel size.
"""
function calc_patch_dimensions(mixer_content::Dict{String,Any})::Tuple{Int64, Int64}
    # Get relevant info from the JSON mixer file
    # NOTE: adding overlap to account for the kernel size,
    # which isn't reflected in the mixer file
    patch_width = mixer_content["patchDimensions"][1]
    patch_height = mixer_content["patchDimensions"][2]
    patch_dimensions = (patch_width, patch_height)
    return patch_dimensions
end


"""
Combines all the patches contained in the file located at infilepath
into a TFRecordDataset. Returns a generator that yields each next example
i.e. patch) in that dataset, parsed into an n_bands x lon x lat numpy array.
"""
function read_tfrecord_file(infilepath::String,
                            dims::Tuple{Int64, Int64},
                            bands::Array{String,1};
                            default_val=DEFAULT_VAL)::OrderedDict{Int64, Array{Float32,3}}
    # empty Dict to hold the arrays
    arrays = Dict()
    # read the TFRecord file
    reader = TFRecordReader(relpath(infilepath))
    # beginning of example printout:
    # Example(Features(Dict{AbstractString,Feature}("constant" => Feature(#undef, FloatList(Float32[
    # loop over the examples
    for (ex_n, example) in enumerate(reader)
        # create an empty Array to hold the data
        arr = Array{Float32}(undef, dims[1], dims[2], 5)
        for (i, band) in enumerate(bands)
            band_vals = example.features.feature[band].float_list.value[1:prod(dims)]
            band_arr = reshape(band_vals, dims) #|> transpose
            # replace the default missing-data val
            replace!(band_arr, default_val=>NaN)
            arr[:,:, i] = band_arr
        end
        arrays[ex_n] = arr 
    end
    # sort the dict by key (to avoid writing files out with patches in diff order)
    return sort(arrays)
end


const (DIMS, CRS, XMIN, YMIN, XRES, YRES,
       PATCHES_PER_ROW, TOT_PATCHES) = get_mixer_info(read_mixer_file(ABS_DATA_DIR))

const FILEPATHS = get_file_paths(ABS_DATA_DIR)


"""
Return an output Dict containing the row, column, and patch numbers,
and outfile paths (as values, each of which are sub-Dicts)
for all files (keys).
"""
function get_row_col_patch_ns_allfiles()::Dict{String,Dict{String,Any}}
    # set the starting row, column, and patch counters
    patch_i = 0
    patch_j = 0
    patch_n = 0

    # get the regex pattern
    patt = Regex("(?<=$(PATT_B4_FILENUM))\\d{5}")

    # assert that list is sorted in ascending numerical order
    # NOTE: if this is not true then my method for tracking the row, col, and
    # patch numbers will break!
    filenums = map(x -> parse(Int32, match(patt, splitdir(x)[2]).match), FILEPATHS)
    filenums_plus1 = filenums .+ 1
    diffs = filenums_plus1 .- filenums
    @assert(all(diffs .== 1), ("Filepaths do not appear to be " * 
                               "in numerical order. \n\t"))

    # make the output Dict
    files_dict = Dict()

    # loop over the input files, get the requisite info (row, col, and patch
    # ns; output file paths), and store in the dict
    for (i, infile) in enumerate(FILEPATHS)
        # create the subdict 
        file_dict = Dict()
        patch_is::Array{Int64,1} = []
        file_dict["patch_is"] = patch_is
        patch_js::Array{Int64,1} = []
        file_dict["patch_js"] = patch_js
        patch_ns::Array{Int64,1} = []
        file_dict["patch_ns"] = patch_ns

        # read the file into a TFRecordReader instance
        dataset = TFRecordReader(infile)

        # loop over the patches in the infile, store the row, col, and patch
        # numbers, then correctly increment them
        # NOTE: an 'example' (TFRecord jargon) is the same as a 'patch' (GEE jargon)
        for example in dataset

            # store nums
            append!(file_dict["patch_is"], patch_i)
            append!(file_dict["patch_js"], patch_j)
            append!(file_dict["patch_ns"], patch_n)

            #increment counters
            patch_n += 1
            if patch_j == PATCHES_PER_ROW - 1
                patch_i += 1
                patch_j = 0
            else
                patch_j += 1
            end
        end

        # add this file's file_dict to the files_dict
        files_dict[infile] = file_dict
    end

    return files_dict
end


const FILES_DICT = get_row_col_patch_ns_allfiles() 

"""
Fuse all output files into a single raster
"""
function fuse_all_files()::Array{Float32, 3}
    # make the fused output array
    ncols::Int64 = PATCHES_PER_ROW * DIMS[2]
    nrows::Int64 = (TOT_PATCHES/PATCHES_PER_ROW) * DIMS[1]
    fused::Array{Float32,3} = ones(nrows, ncols, 5) * NaN

    # loop over output files, reading each one in and putting its data into the fused array
    for (filename, fileinfo) in FILES_DICT
        patches = read_tfrecord_file(filename, DIMS, BANDS)
        patch_nums = sort(collect(keys(patches)))

        # loop over patches in the file
        for patch_num in patch_nums
            # get corresponding indices in the fused array
            i_start = fileinfo["patch_is"][patch_num]*DIMS[1] + 1
            i_stop = i_start + DIMS[1] - 1
            j_start = fileinfo["patch_js"][patch_num]*DIMS[2] + 1
            j_stop = j_start + DIMS[2] - 1

            # add data to the fused array
            fused[i_start:i_stop, j_start:j_stop, :] = patches[patch_num]
        end

    end
       
    return fused

end

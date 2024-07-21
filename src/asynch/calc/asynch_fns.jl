"""
Reads in the data from all the TFRecord files pertaining to a mixerfile.
Calculates asynchrony for that data and writes out to a matching set of files.

NOTES:
    - TFRecord terminology:
        Each TFRecord file is composed of multiple 'examples'
        (GEE equiv: 'patches'), each of which is composed of multiple
        'features' (GEE: 'bands'), each of which has both
        a 'key' (GEE: 'band name', e.g. 'sin_1') and a value
        (GEE's actual raster of image data; in TFRecord they are stored
        and reading in as 'float_lists' composed of numerous 'value: #######'
        pairs, where ##### is the actual numerical value in a cell).

        In other words, the TFRecord format looks like:

            file_n.tfrecords:

                example 1:

                    feature {
                        key: "constant"
                        value {
                            float_list {
                                value: -9999.0
                                value: -9999.0
                                value: -9999.0
                                .
                                .
                                .
                                }
                            }
                        key: "sin_1"
                        ...
                        key: "sin_2"
                        ...
                        ...
                        }

                example 2:

                    feature {
                        ...
                        }
                ...

                example last:

                    feature {
                        ...
                        }


"""

using OrderedCollections
using NearestNeighbors
using StaticArrays
using Distributed
using Statistics
using StatsBase
using Distances
#using GeoArrays
using TFRecord
using ArchGDAL
using Printf
#using Images
using Colors
#using Plots
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
const VERBOSE = false
"""
Whether to time the asynchrony calculation
(Time will be displayed per pixel)
"""
const TIMEIT = false

# either get path on my laptop...
if splitpath(pwd())[2] == "home"
    const BASE_DATA_DIR = "/media/deth/SLAB/diss/3-phn/GEE_outputs/"
# ... or else get path on Savio
else
    const BASE_DATA_DIR = "/global/scratch/users/drewhart/seasonality/quick_GEE_outputs/"
end

"""
indicator of whether or not TFRecordReader is defined
"""
# code to handle differing versions of Julia, and hence of TFRecord,
# on Savio and on my laptop
const HAS_TFRECORDREADER = @isdefined TFRecordReader

"""
pattern that occurs just before the file number
"""
const PATT_AFT_FILENUM = "\\.tfrecord\$"

"""
kernel size used by GEE to output the TFRecord files
"""
# kernel size
const KERNEL_SIZE = 288
 
"""
half the kernel width (to use in margin-trimming)
"""
const HKW = floor(Int, KERNEL_SIZE/2)

"""
whether or not to trim the half-kernel-width margin before output
"""
const TRIM_MARGIN = true

"""
default missing-data val
"""
const DEFAULT_VAL = -9999.0

"""
names of the bands saved in the input TFRecord files
NOTE: ignores the last band, which is the R^2 values
"""
const INBANDS = ["constant", "sin_1", "cos_1", "sin_2", "cos_2"]

"""
names of the bands to save in the output TFRecord files
"""
const OUTBANDS = ["asynch", "asynch_r2", "asynch_pval", "asynch_n"]

"""
approximate radius of the earth (for the Haversine distance-based,
second-pass neighbor filtering)
NOTE: 6371008 m is a good approximation of any of the three common methods for
      representing the Earth by its 'mean' radius
"""
const EARTH_RAD = 6_371_008

"""
minimum number of neighbors needed for a pixel's asynchrony to be calculated
and for its result to be included in the output dataset
"""
const MIN_NUM_NEIGHS = 30

"""
whether or not to fit interecepts in the asynchrony-calculating regressions
"""
const FIT_INTERCEPT = true

"""
α of the regression slopes
(i.e., max slope P-value to consider significant)
"""
const ALPHA = 0.01


#-----------------
# define functions
#-----------------
#
"""
Uses the data-directory path provided to get and return
lists of the input and output files' paths.
"""
function get_infile_outfile_paths(data_dir::String,
                                  neigh_rad::Int64)::Tuple{Array{String,1}, Array{String,1}}
    # set up IO paths
    infilepaths = glob("*tfrecord", data_dir) 
    # order the infilepaths
    sort!(infilepaths)
    # exclude any previously generated output files
    infilepaths = [fp for fp in infilepaths if !occursin("-OUT-", fp)]
    # get the corresponding outfilepaths
    outfilepaths = [replace(fp, r"-(?=\d{5})" => string(neigh_rad) * "mrad-OUT-") for fp in infilepaths]
    return infilepaths, outfilepaths
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
    dims = calc_patch_dimensions(mixer_content, KERNEL_SIZE)
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
    xmin = affine[3] - (((KERNEL_SIZE/2) - 0.5)* xres)
    ymin = affine[6] - (((KERNEL_SIZE/2) - 0.5)* yres)
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
function calc_patch_dimensions(mixer_content::Dict{String,Any},
                               kernel_size::Int64)::Tuple{Int64, Int64}
    # Get relevant info from the JSON mixer file
    # NOTE: adding overlap to account for the kernel size,
    # which isn't reflected in the mixer file
    patch_width = mixer_content["patchDimensions"][1] + kernel_size
    patch_height = mixer_content["patchDimensions"][2] + kernel_size
    patch_dimensions = (patch_width, patch_height)
    return patch_dimensions
end


"""
Combines all the patches contained in the file located at infilepath
into a TFRecordDataset. Returns an integer-keyed dict of
examples (i.e. patches) each one being an n_bands x lon x lat array.
"""
function read_tfrecord_file(infilepath::String,
                            dims::Tuple{Int64, Int64},
                            bands::Array{String,1};
                            default_val=DEFAULT_VAL)::OrderedDict{Int64, Array{Float32,3}}
    # empty Dict to hold the arrays
    arrays = Dict()
    # read the TFRecord file
    # beginning of example printout:
    # Example(Features(Dict{AbstractString,Feature}("constant" => Feature(#undef, FloatList(Float32[
    # loop over the examples
    if HAS_TFRECORDREADER
        for (ex_n, example) in enumerate(TFRecordReader(infilepath))
            # create an empty Array to hold the data
            arr = Array{Float32}(undef, dims[1], dims[2], 5)
            for (i, band) in enumerate(bands)
                band_vals = example.features.feature[band].float_list.value[1:prod(dims)]
                # NOTE: need to transpose array because tfrecord store in row-major order
                #       but julia's reshape fills in column-major order
                band_arr = transpose(reshape(band_vals, dims))
                # replace the default missing-data val
                replace!(band_arr, default_val=>NaN)
                arr[:,:, i] = band_arr
            end
            arrays[ex_n] = arr 
        end
    else
        for (ex_n, example) in enumerate(TFRecord.read(infilepath))
            # create an empty Array to hold the data
            arr = Array{Float32}(undef, dims[1], dims[2], 5)
            for (i, band) in enumerate(bands)
                band_vals = example.features.feature[band].float_list.value[1:prod(dims)]
                # NOTE: need to transpose array because tfrecord store in row-major order
                #       but julia's reshape fills in column-major order
                band_arr = transpose(reshape(band_vals, dims))
                # replace the default missing-data val
                replace!(band_arr, default_val=>NaN)
                arr[:,:, i] = band_arr
            end
            arrays[ex_n] = arr 
        end
    end
    # sort the dict by key (to avoid writing files out with patches in diff order)
    return sort(arrays)
end


"""
Writes the output patches (i.e. rasters) to disk as a TFRecord file,
using the given filepath.
"""
function write_tfrecord_file(patches::OrderedDict{Int64, Array{Float32,3}},
			     filepath::String,
                             bands::Array{String, 1};
                             default_val=DEFAULT_VAL)::Nothing
    #instantiate the TFRecordWriter
    writer = TFRecord.TFRecordWriter(filepath)
    # NOTE: sort the patches once more, to ensure they're in the right order
    for (patch_i, patch) in sort(patches)
        # create a Dict of the outbands and their arrays
        outdict = Dict()
        for (band_i, band) in enumerate(bands)
            # set all missing back to the missing-data default val,
            # then recast as Float32 and cast as a vector
            # NOTE: need to undo the transpose we did when reading the TFRecord file in,
            #       so that julia's column-major behavior returns data in the row-major
            #       order used in the original TFRecord file
            outpatch = vec(Float32.(replace(transpose(patch[:, :, band_i]), NaN=>default_val)))
            outdict[band] = outpatch
        end
        # serialize to a TF example
        write(writer, outdict)
    end
    close(writer)
end


"""
Preps a single patch for writing into a TFRecord file
"""
function prep_single_patch_for_tfrecord_file(patch::Array{Float32,3},
                                             bands::Array{String, 1};
                                             default_val=DEFAULT_VAL)::Dict{String,
                                                                           Vector{Float32}}
    patchdict = Dict()
    for (band_i, band) in enumerate(bands)
        # NOTE: need to undo the transpose we did when reading the TFRecord file in,
        #       so that julia's column-major behavior returns data in the row-major
        #       order used in the original TFRecord file
        outpatch = vec(Float32.(replace(transpose(patch[:, :, band_i]), NaN=>default_val)))
        patchdict[band] = outpatch
    end
    return patchdict
end


"""
Writes the output patches (i.e. rasters) to disk as a TFRecord file,
using the given filepath.
"""
function write_tfrecord_file_new(patches::OrderedDict{Int64, Array{Float32,3}},
                             filepath::String,
                             bands::Array{String, 1})::Nothing
    TFRecord.write(filepath, (prep_single_patch_for_tfrecord_file(patch_item[2], bands) for patch_item in sort(patches)))
end


"""
Takes the overall x and y min values of the TFRecord dataset,
the x and y resolutions, the patch dimensions, and the column and row
indices of the current patch. Returns x and y grids of the
cell centers of the current patch.
"""
function get_patch_lons_lats(xmin::Float64, ymin::Float64, xres::Float64, yres::Float64,
                             dims::Tuple{Int64,Int64}, hkw::Int64,
                             patch_j::Int64, patch_i::Int64)::Tuple{Array{Float64,2},Array{Float64,2}}
    # calculate the x and y mins of the current patch
    kern_xdim, kern_ydim = dims
    real_xdim, real_ydim = [kern_xdim kern_ydim] .- (KERNEL_SIZE)
    patch_xmin = xmin + (patch_j * real_xdim * xres)
    patch_ymin = ymin + (patch_i * real_ydim * yres)
    # NOTE: -1 gets us from the center of the first cell to the center of the last
    patch_xmax = patch_xmin + (xres * (kern_xdim-1))
    patch_ymax = patch_ymin + (yres * (kern_ydim-1))

    # get lists of xs and ys of all the current patch's pixels
    
    xs = LinRange(patch_xmin, patch_xmax, kern_xdim)
    ys = LinRange(patch_ymin, patch_ymax, kern_ydim)
    # get the meshgrids of those coordinates
    gridx = xs' .* ones(length(xs))
    gridy = ys .* ones(length(ys))'
    # check that y values are identical across rows and x values are identical down columns
    @assert(unique([length(unique(gridy[i,:])) == 1 for i in 1:size(gridy)[1]]) == [1,],
            "y values are not identical across rows in gridy!")
    @assert(unique([length(unique(gridx[:,j])) == 1 for j in 1:size(gridx)[2]]) == [1,],
            "x values are not identical down columns in gridx!")
    return gridx, gridy
end


"""
Makes and returns the regression's design matrix, a 365 x 5 array
in which the columns contain, in order:
    - 1s (for the constant);
    - sin and cos of annual-harmonic days of the year
      (i.e. days are expressed in radians from 0 to 2pi);
    - sin and cos of the semiannual-harmonic days of the year
      (i.e. days are expressed in radians from 0 to 4pi).
"""
function make_design_matrix()::Array{Float64,2}
    # get 1 year of daily values, expressed in radians, 1 rotation/yr
    annual_radian_days = LinRange(0, 2pi, 366)[1:365]
    # get 1 year of daily values, expressed in radians, 2 rotations/yr
    semiannual_radian_days = @. LinRange(0, 4pi, 366)[1:365] % (2pi)
    # get the harmonic values of those
    sin1 = sin.(annual_radian_days)
    cos1 = cos.(annual_radian_days)
    sin2 = sin.(semiannual_radian_days)
    cos2 = cos.(semiannual_radian_days)
    # add a vector of 1s for the constant term, then recast as a 365 x 5 array,
    # to use as the covariate values in the regression
    design_mat = [ones(size(sin1)) sin1 cos1 sin2 cos2]
    return design_mat
end


"""
Calculates the time series at pixel i,j, using the coefficients for the
constant and the sin and cosine terms of the annual and semiannual
harmonic components. Returns the time series as an array.
"""
function calc_time_series(patch::Array{Float32,3}, i::Int64, j::Int64,
                          design_mat::Array{Float64,2})::Array{Float64,1}
    # multiply the pixel's set of coefficients by the design mat, then sum
    # all the regression terms
    # NOTE: the coeffs are array of shape (5,);
    #       the design matrix is array of shape (365, 5);
    #       the multiplication operates row by row, and elementwise on each
    #       row, returning a new (365, 5) array that, when summed across the
    #       columns (i.e. the regression's linear terms) gives a (365,) vector
    #       of fitted daily values for the pixel
    ts = vec(sum(patch[i, j, :]' .* design_mat, dims=2))
    return ts
end


"""
Calculates the simple linear regression of y on x,
where variables are column names in DataFrame df.

Returns the intercept, slope, and R^2 in a Dict.
"""
function run_linear_regression(x::Array{Float64,1}, y::Array{Float64,1};
                               fit_intercept=true)::Dict{String, Float64}
    # make the design matrix
    if fit_intercept
        X = hcat(ones(length(x)), x)
    else
        X = reshape(x, length(x), 1)
    end

    # fit a model
    try
        mod = fit(LinearModel, X, y)

        # get the R2 and coeff values
        R2 = r2(mod)
        if fit_intercept
            int, slope = coef(mod)
            Pval = coeftable(mod).cols[4][2]
        else
            int = NaN
            slope = coef(mod)[1]
            Pval = coeftable(mod).cols[4][1]
        end
        output = Dict("slope" => slope,
                      "intercept" => int,
                      "R2" => R2,
                      "Pval" => Pval)
        return output
    catch error
        R2 = NaN
        int = NaN
        slope = NaN
        Pval = NaN
        output = Dict("slope" => slope,
                      "intercept" => int,
                      "R2" => R2,
                      "Pval" => Pval)
        @warn "ERROR THROWN!"
        return output
    end
end


"""
Standardizes a 1d Array. Returns the standardized array of identical shape.
"""
standardize(a1::Array{Float64,1})::Array{Float64,1} = (a1.-mean(a1))/std(a1)


"""
Uses the dataset's mixer info to create a lookup dict
containing the number of neighors within the neigh_rad-radius
neighborhood for each latitudinal band of cells.
"""
function make_nneighs_lookup_dict(mix_info::Dict{String, Any},
                                  hkw::Int64,
                                  neigh_rad::Int64)::Dict{Tuple{Int64, Int64}, Int64}

    # get min y value (i.e. northmost)
    # and yres (which will be negative, hence progressing southward)
    ymin = mix_info["projection"]["affine"]["doubleMatrix"][6]
    yres = mix_info["projection"]["affine"]["doubleMatrix"][5]

    # get number of patch rows
    # and total number of lat cells
    n_patch_rows = mix_info["totalPatches"]/mix_info["patchesPerRow"]
    @assert(n_patch_rows%1 == 0, "n_patch_rows is not an integer!")
    n_patch_rows = floor(Int, n_patch_rows)

       
    patch_ydim = mix_info["patchDimensions"][2]
    tot_lat_cells = n_patch_rows * patch_ydim
    # NOTE: should already be an integer, but can't hurt to double-check
    @assert(tot_lat_cells%1 == 0, "tot_lat_cells is not an integer!")
    tot_lat_cells = floor(Int, tot_lat_cells)

    # get all lat values
    # NOTE: add half res so that we're dealing with cell centers, not cell corners
    lats = LinRange(ymin + 0.5*yres,
                    ymin + 0.5*yres + (yres * (tot_lat_cells - 1)),
                    tot_lat_cells)

    # get the pairs of patch row-number and cell row-number
    # that correspond to each cell-row in the final output dataset
    # NOTE: have to be 0-indexed in order to match up with
    #       the patch_is values being returned by get_row_col_patch_ns_allfiles
    patch_is = repeat([0:n_patch_rows-1;], inner=patch_ydim)
    # NOTE: add hkw so that the is are expressed as they will appear in the overlapping
    #       patch arrays, before we trim the margins!
    row_is = repeat([1+hkw:patch_ydim+hkw;], n_patch_rows)

    # pair up actual cell-center lats with their patch and cell i values
    # (as a Dict, with Tuple(patch_i, cell_i) as keys and lat as values)
    cells_lats_dict = Dict(zip(zip(patch_is, row_is), lats))

    # create an empty Dict to fill up with our output
    # (with Tuple(patch_i, cell_i) as keys and nneighs as values)
    nneighs_lookup_dict = Dict()
 
    # get xres and patch xdim (for use creating an arbitrary
    # range of potential longitude values in the loop below)
    xres = mix_info["projection"]["affine"]["doubleMatrix"][1]

    # for each item in the cells_lats_dict
    for (cell_info, lat) in cells_lats_dict
        # NOTE: add 1 to the total number of cells so that the focal cell,
        # whose center coords are (0, lat), will be the center cell
        # in the potential-neighbor box we create
        potent_lons = LinRange(0-(2*hkw*xres), 0+(2*hkw*xres), 4*hkw+1)
        potent_lats = LinRange(lat-(2*hkw*yres), lat+(2*hkw*yres), 4*hkw+1)
        vec_potent_lons = vec(potent_lons' .* ones(length(potent_lons)))
        vec_potent_lats = vec(potent_lats .* ones(length(potent_lats))')
        # build a tree containing excess potential coords
        tree = BallTree([vec_potent_lons vec_potent_lats]', Haversine(EARTH_RAD))
        # query the tree and record number of neighs
        nneighs = length(inrange(tree, [0, lat], neigh_rad))
        # store the number of neighs in the cells_nneighs_dict
        nneighs_lookup_dict[cell_info] = nneighs
    end
    
    return nneighs_lookup_dict

end


"""
Finds all pixel-center coordinates within neigh_dist (in meters) of
the input lat,lon coordinates.

To do this, a BallTree (constructed to use the Haversine distance metric)
is first queried to find all neighbors within the neighborhood radius (in meters)
(rather than my former approach, which was based on the
minimum distance of a degree longitude at my highest latitude, i.e. 60° N/S)

Then a Haversine distance calculator is used to calculate the actual distance,
in meters, of each pixel's centerpoint from the focal pixel's centerpoint.
NOTE: Could use the knn function called on the BallTree instead, but quick profiling at
the REPL shows that to be ~175x slower (0.004068 sec vs 0.000023 sec, for 874 neighs,
using the @time macro)

Neighbor pixels' (lat,lon coordinates) => (distances, (array indices)) are returned
as key => value pairs in a dict.
"""
function get_neighbors_info(i::Int64, j::Int64,
                            vec_is::Array{Int64,1}, vec_js::Array{Int64,1},
                            foc_y::Float64, foc_x::Float64,
                            vec_ys::Array{Float64}, vec_xs::Array{Float64},
                            yres::Float64, xres::Float64,
                            patch_dims::Tuple{Int64,Int64},
                            patch_i::Int64, nneighs_lookup_dict::Dict{Tuple{Int64, Int64}, Int64},
                            hav_fn::Function,
                            tree::BallTree{SArray{Tuple{2},Float64,1,2},2,Float64,Haversine{Int64}})::Dict{Tuple{Float64, Float64}, Tuple{Float64, Tuple{Int64, Int64}}}
    # get the tree's data-row numbers for all neighbors
    # within the neighborhood (by grabbing the k nearest neighbors, where k comes
    # from the nneighs_lookup_dict built at the outset with neigh_rad factored in)
    neighs = knn(tree, [foc_x, foc_y], nneighs_lookup_dict[Tuple((patch_i, i))])[1]
    # use the row numbers to subset the pixels' centerpoint coords and
    # their index coordinate pairs
    neigh_ys = vec_ys[neighs]
    neigh_xs = vec_xs[neighs]
    neigh_is = vec_is[neighs]
    neigh_js = vec_js[neighs]
    # make sure num entities in neigh_coords == num index-pairs in neigh_indices
    @assert(length(neigh_xs) == length(neigh_js), ("Incorrect numbers of " *
                                                   "valid x coords and/or " *
                                                   "j indices returned!"))
    @assert(length(neigh_ys) == length(neigh_is), ("Incorrect numbers of " *
                                                   "valid y coords and/or " *
                                                   "i indices returned!"))

    # calculate the haversine distances
    # NOTE: the Haversine formula is less accurate than the Vincenty,
    #       but also less computationally intensive, and thus faster;
    #       using it for now, since compute time is my biggest concern by far,
    #       and since the error is quite small (appears to be on the order of 
    #       a few hundred meters for points up to ~150km apart, so this should
    #       not be a problem for the ~5 km-res SIF data, but could later
    #       present an issue if I decide to crunch the ~500 m-res NIRv data...)
    #
    neigh_dists = hav_fn.(zip(neigh_xs, neigh_ys))

    # combine into a dict
    out = Dict(zip(zip(neigh_xs, neigh_ys),
                   # NOTE: GEOGRAPHIC COORDS ARE x,y ORDER, BUT CELL COORDS ARE i,j!
                   zip(neigh_dists, zip(neigh_is, neigh_js))))
    return out
end


"""
Makes vectors of the i and j indices of the patch's pixels.
The structure and order of these vectors makes them
subsettable by the output from a neighbor search on a BallTree.
"""
function get_is_js(dims::Tuple{Int64, Int64})::Tuple{Array{Int64,1}, Array{Int64,1}}
    is = vec(repeat([1:dims[1];], dims[2]))
    js = vec(repeat([1:dims[2];], inner=dims[1]))
    return is, js
end


"""
Read the file's data from a TFRecord file, as a set of Arrays,
one per patch.

Return both that data and a set of identically structured output Arrays.
"""
function get_inpatches_outpatches(infilepath::String, inbands::Array{String,1},
         dims::Tuple{Int64, Int64})::Tuple{OrderedDict{Int64, Array{Float32,3}},
                                           OrderedDict{Int64, Array{Float32,3}}}

    # read the file's data in as a set of examples (i.e. patches)
    inpatches = read_tfrecord_file(infilepath, dims, inbands)
    # create a data container for all the asynch patches
    outpatches = Dict()
    for (i, patch) in inpatches
        # create the output arrays, to be filled in (start as all NaNs),
        # but specify the same sizes on dims 1 and 2 and new output size on dim 3
        # (dim 3 is size 4 because it holds layers for the asynch metric, its R^2, its P-value, and its n)
        outpatches[i] = similar(patch, size(patch)[1], size(patch)[2], 4) * NaN
    end
    return (inpatches, sort(outpatches))
end


"""
Return an output Dict containing the row, column, and patch numbers,
and outfile paths (as values, each of which are sub-Dicts)
for all files (keys).
"""
function get_row_col_patch_ns_allfiles(data_dir::String,
                                       patt_aft_filenum::String,
                                       neigh_rad::Int64)::Dict{String,Dict{String,Any}}
    # set the starting row, column, and patch counters
    patch_i = 0
    patch_j = 0
    patch_n = 0

    # get the mixer file info
    mix = read_mixer_file(data_dir)
    (dims, crs, xmin, ymin, xres, yres,
     patches_per_row, tot_patches) = get_mixer_info(mix)

    # get all the input and output file paths
    infilepaths, outfilepaths = get_infile_outfile_paths(data_dir, neigh_rad)

    # get the regex pattern
    patt = Regex("\\d{5}?(?=$patt_aft_filenum)")

    # assert that both lists are sorted in ascending numerical order
    # NOTE: if this is not true then my method for tracking the row, col, and
    # patch numbers will break!
    for filepaths in [infilepaths, outfilepaths]
        filenums = map(x -> parse(Int32, match(patt, splitdir(x)[2]).match),
                       filepaths)
        filenums_plus1 = filenums .+ 1
        diffs = filenums_plus1 .- filenums
        @assert(all(diffs .== 1), ("Filepaths do not appear to be " * 
                                   "in numerical order. \n\t"))
    end

    # make the output Dict
    files_dict = Dict()

    # loop over the input files, get the requisite info (row, col, and patch
    # ns; output file paths), and store in the dict
    for (infile_i, infile) in enumerate(infilepaths)
        # create the subdict 
        file_dict = Dict()
        # stash the outfile path
        file_dict["outfilepath"] = outfilepaths[infile_i]
        patch_is::Array{Int64,1} = []
        file_dict["patch_is"] = patch_is
        patch_js::Array{Int64,1} = []
        file_dict["patch_js"] = patch_js
        patch_ns::Array{Int64,1} = []
        file_dict["patch_ns"] = patch_ns

        # read the file into a TFRecordReader instance
        if HAS_TFRECORDREADER
            dataset = TFRecordReader(infile)
        else
            dataset = TFRecord.read(infile)
        end
        #opened_file = TFRecord.read(infile)

        # loop over the patches in the infile, store the row, col, and patch
        # numbers, then correctly increment them
        # NOTE: an 'example' (TFRecord jargon) is the same as a 'patch' (GEE jargon)
        for example in dataset
        #for example in opened_file

            # store nums
            append!(file_dict["patch_is"], patch_i)
            append!(file_dict["patch_js"], patch_j)
            append!(file_dict["patch_ns"], patch_n)

            #increment counters
            patch_n += 1
            if patch_j == patches_per_row - 1
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


"""
Calculates and sets the asynchrony value for a single pixel
located at position i,j in the patch.
"""
function calc_asynch_one_pixel!(i::Int64, j::Int64,
                               vec_is::Array{Int64,1}, vec_js::Array{Int64,1},
                               foc_y::Float64, foc_x::Float64,
                               vec_ys::Array{Float64,1}, vec_xs::Array{Float64,1},
                               patch::Array{Float32,3}, outpatch::Array{Float32,3},
                               patch_n::Int64, yres::Float64, xres::Float64,
                               dims::Tuple{Int64, Int64},
                               design_mat::Array{Float64,2},
                               patch_i::Int64, nneighs_lookup_dict::Dict{Tuple{Int64, Int64}, Int64},
                               tree::BallTree{SArray{Tuple{2},Float64,1,2},2,Float64,Haversine{Int64}};
                               timeit=true, verbose=true)::Nothing
    if verbose && timeit
        # get start time
        start = now()
    end

    if verbose
        println("\nPROCESSING PATCH $patch_n: PIXEL ($i, $j)...")
        println("\tLAT: $foc_y; LON: $foc_x")
    end

    # create lists of ts dist and geo dist values
    ts_dists::Array{Float64,1} = []
    geo_dists::Array{Float64,1} = []

    # calculate the focal pixel's time series (and its standardized time series)
    ts_foc = calc_time_series(patch, i, j, design_mat)
    stand_ts_foc = standardize(ts_foc)

    # make the Haversine instance and function for this cell
    hav_inst = Haversine(EARTH_RAD)
    hav_fn = function(coord_pair)
        return hav_inst((foc_x, foc_y), coord_pair)
    end

    # get the coords, dists, and array-indices
    # of all of the focal pixel's neighbors
    coords_dists_inds = get_neighbors_info(i, j, vec_is, vec_js,
                                           foc_y, foc_x, vec_ys, vec_xs,
                                           yres, xres, dims,
                                           patch_i, nneighs_lookup_dict,
                                           hav_fn, tree)

    if verbose
    	println("\t\t($foc_x, $foc_y)\n")
    	println("\t\t$(length(coords_dists_inds)) neighbors\n")
    end

    # loop over neighbors
    for (neigh_coords, neigh_info) in coords_dists_inds
        # unpack the neighbor info
        neigh_dist, neigh_inds = neigh_info
        ni, nj = neigh_inds

        # get the neighbor's time series
        ts_neigh = calc_time_series(patch, ni, nj, design_mat)
        stand_ts_neigh = standardize(ts_neigh)

        # drop this pixel if it returns NAs
        if any(isnan.(view(patch, ni,nj,:)))
           # do nothing 
           continue
        else
            try
                # calculate the Euclidean distance
                # between the 2 standardized time series
                ts_dist = euclidean(stand_ts_foc, stand_ts_neigh)
         
                # append the metrics, if the neighbor-distance
                # calculation was successful
                append!(geo_dists, neigh_dist)
                append!(ts_dists, ts_dist)
            catch error
                @warn ("ERROR THROWN WHILE CALCULATING NEIGHBOR-DISTANCE METRICS:" *
                         "\n\tpixel: $i, $j\n\tneighbor: $ni, $nj\n")
            end
        end
    end

    # calc asynchrony for the cell, but only if we have enough valid neighbors
    num_neighs = length(geo_dists)
    if verbose
        println("$num_neighs NEIGHBORS")
    end
    if num_neighs ≥ MIN_NUM_NEIGHS
        # get the slope of the overall regression of Euclidean ts dists on geo dist
        # NOTE: setting fit_intercept to false fits ts_dist to 0 at geo_dist=0
        res_ts = run_linear_regression(geo_dists, ts_dists, fit_intercept=FIT_INTERCEPT)
        # extract both R2 and Pval into vars
        R2_ts = res_ts["R2"]
        Pval_ts = res_ts["Pval"]
        # extract slope (i.e., asynchrony metric),
        # and use Pval to determine if it is indistinguishable from zero
        asynch_ts = res_ts["slope"] * Int(Pval_ts < ALPHA)
        # extract sample size (i.e. number of neighbors) for this focal pixel,
        # checking that there was no issue with uneven numbers of results
        @assert(length(geo_dists) == length(ts_dists)), ("Lengths of " *
                                                                       "geo_dists, " *
                                                                       "ts_dists, " *
                                                                       "and R2s are " *
                                                                       "not equal!")

        # update the output patch with the new values
        outpatch[i, j, :] = [asynch_ts, R2_ts, Pval_ts, convert(Float32, num_neighs)]

    # if we don't have at least the minimum number of neighbors
    # then just add NaNs to the outpatch, except for the neighbor count
    else
        outpatch[i, j, :] = [NaN, NaN, NaN, convert(Float32, num_neighs)]
        if verbose
            println("TOO FEW NEIGHBORS!")
	end
    end

    # calculate runtime and print, if necessary
    if verbose && timeit
        # get finish time, print out total time
        stop = now()
        diff = stop-start
        println("\n\truntime: $diff")
        println("==================================")
    end
end


"""
Read the patches from the input file, and calculate and store
the asynch metrics for each input patch (inpatches)
in the corresponding output patch (outpatches)
"""
function calc_asynch(inpatches::OrderedDict{Int64, Array{Float32,3}},
                     outpatches::OrderedDict{Int64, Array{Float32,3}},
                     patch_is::Array{Int64,1}, patch_js::Array{Int64,1}, patch_ns::Array{Int64,1},
                     vec_is::Array{Int64,1}, vec_js::Array{Int64,1},
                     cart_inds::CartesianIndices{2,Tuple{Base.OneTo{Int64},Base.OneTo{Int64}}},
                     xmin::Float64, ymin::Float64, xres::Float64, yres::Float64,
                     dims::Tuple{Int64,Int64}, design_mat::Array{Float64,2},
                     nneighs_lookup_dict::Dict{Tuple{Int64, Int64}, Int64},
                     kernel_size::Int64, hkw::Int64;
                     trim_margin=false, verbose=true, timeit=true)::OrderedDict{Int64, Array{Float32,3}}

    # loop over the patches from the current file
    # NOTE: each patch is of shape (n_bands, lat, lon)
    for (patch_idx, inpatch) in inpatches
        # grab important variables
        outpatch = outpatches[patch_idx]
        patch_i = patch_is[patch_idx]
        patch_j = patch_js[patch_idx]
        patch_n = patch_ns[patch_idx]

        # get the lons and lats of the pixels in the current example's patch
        # as well as an array containing columns for each of the coord pairs of the pixels
        xs, ys = get_patch_lons_lats(xmin, ymin, xres, yres, dims, hkw, patch_j, patch_i)
        vec_ys = vec(ys)
        vec_xs = vec(xs)

        # get the BallTree for these lons and lats
        # make a KDTree out of all the coordinates in a patch,
        # to be used in neighbor-searching
        tree = BallTree([vec_xs vec_ys]', Haversine(EARTH_RAD))

        #----------------
        # run calculation
        #----------------
  
        # data structs to make sure that each foc_x and foc_y are hit equal number of times
        # (i.e. that there's not procession in these values as we progress across cols)
        foc_pts = []

        # loop over pixels (excluding those in the kernel's margin,
        # since there's no use wasting time calculating for those)
        for ind in eachindex(inpatch[:, :, 1])
            i,j  = Tuple(cart_inds[ind])
            if (hkw < i <= (size(inpatch)[1]-hkw)) && (hkw < j <= (size(inpatch)[2]-hkw))

                # leave the asynch output val as a NaN if the coeffs for this
                # pixel contain NaNs
                if any(isnan.(inpatch[i, j, :]))
                    if verbose
                        println("\nPROCESSING PATCH $patch_n: PIXEL ($i, $j)...")
                        lat = ys[i, j]
                        lon = xs[i, j]
                        println("\tLAT: $lat; LON: $lon")
                        println("\n\truntime: ---")
                        println("==================================")
                    end
                else
                    foc_y, foc_x = (ys[i,j], xs[i,j])
                    push!(foc_pts, (foc_x, foc_y))
                    calc_asynch_one_pixel!(i, j, vec_is, vec_js,
                                          foc_y, foc_x, vec_ys, vec_xs,
                                          inpatch, outpatch, patch_n,
                                          yres, xres, dims,
                                          design_mat,
                                          patch_i, nneighs_lookup_dict,
                                          tree, verbose=verbose, timeit=timeit)
                end
            else
                if verbose
                    println("        skipping $i, $j")
                end
            end
        end

        # make sure that all foc_x and foc_y values were visited equally often
        foc_pt_cts = countmap(foc_pts)
        # NOTE: make <=1 so that length 0 (when no non-nan focal pts exist)
        #       doesn't break the code
        @assert(length(unique(values(foc_pt_cts))) <= 1, "not all foc_pts occurred equally!\n\n$foc_pts\n\n")
        #@assert(length(foc_pt_cts) == unique(values(foc_pt_cts))[1], "number of foc_pts values used not equal to number of time each point was used!")

    end

    # trim the half-kernel-width margin, if requested
    if trim_margin
        if verbose
            println("\n\nTrimming $hkw-cell margin around each patch.\n\n")
        end
        # subset
        trimmed_outpatches = Dict()
        for (i, op) in outpatches
            trimmed_outpatches[i] = op[hkw+1:end-hkw, hkw+1:end-hkw, :]
        end
        return trimmed_outpatches
    else
        return deepcopy(outpatches)
    end
end



#----------------------------------------------------
# define the main function, to be mapped over workers
#----------------------------------------------------

"""
Main function, to be mapped over multiple workers.

Takes the file_info for an input file, calculates
asynchrony for that file, writes it to an output TFRecord file,
returns nothing.

The argument 'file_info' must be a Dict item object of the form:
  Dict(infilepath =>
        {"outfilepath": outfilepath,
         "patch_is" => [patch_i_1, patch_i_2, ..., patch_i_I],
         "patch_js" => [patch_j_1, patch_j_2, ..., patch_j_J],
         "patch_ns" => [patch_n_1, patch_n_2, ..., patch_n_N]
        }
      )
"""
function main_fn(file_info::Tuple{String,Int64,String,Dict{String,Any}};
                 verbose=VERBOSE, timeit=TIMEIT, trim_margin=TRIM_MARGIN)::Nothing

    # split input info into variables
    data_dir, neigh_rad, infilepath, file_dict = file_info

    # read in information from the JSON mixer file.
    mix = read_mixer_file(data_dir)
    (dims, crs, xmin, ymin, xres, yres,
     patches_per_row, tot_patches) = get_mixer_info(mix)
    
    # make dict containing (patch_i, row_i) as keys and
    # number of neighbor cells as values;
    # will be used to subset the correct number (k) of nearest neighbors
    # from the Haversine BallTree for each focal cell
    nneighs_lookup_dict = make_nneighs_lookup_dict(mix, HKW, neigh_rad)
    
    # get vectors of pixel indices in the i and j array dims
    # (which will be subsetted by BallTree neighbor-query output)
    vec_is, vec_js = get_is_js(dims)
    
    # a CartesianIndices object to be used for iteration over
    # a patch's pixels in calc_asynchrony();
    # tried this approach while reading about Julia Array types,
    # to help learn better how they are put together and work,
    # and while it appears to only offer a very slight speed-up,
    # it surely can't hurt
    cart_inds = CartesianIndices(dims)

    # make the design matrix from the harmonic regression,
    # to be paired with the coefficients in order to
    # calculate the 365-day fitted time series
    design_mat = make_design_matrix()

    outfilepath = file_dict["outfilepath"]
    patch_is::Array{Int64,1} = file_dict["patch_is"]
    patch_js::Array{Int64,1} = file_dict["patch_js"]
    patch_ns::Array{Int64,1} = file_dict["patch_ns"]
    if verbose
        @info "\nWorker $(myid()) processing file $infilepath\n"
    end

    # read the data in, and set up the output patches' data structure
    # NOTE: by setting this up outside the calc_asynch function, 
    #       I make it so that I can run the script for a short amount of
    #       time, then retain the partial result
    #       for interactive introspection
    inpatches, outpatches = get_inpatches_outpatches(infilepath, INBANDS, dims)

    if verbose
        println("RUNNING ASYNCH CALC FOR FILE: $infilepath")
        for (patch_i, patch_j, patch_n) in zip(patch_is, patch_js, patch_ns)
            println("\tPATCH: $patch_n (patch row: $patch_i, patch col: $patch_j)")
        end
    end

    # run the asynchrony calculation
    outpatches = calc_asynch(inpatches, outpatches,
                             patch_is, patch_js, patch_ns, vec_is, vec_js, cart_inds,
                             xmin, ymin, xres, yres, dims, design_mat, nneighs_lookup_dict, KERNEL_SIZE, HKW,
                             trim_margin=trim_margin, verbose=verbose, timeit=timeit)

    # write out the asynch data
    if HAS_TFRECORDREADER
        write_tfrecord_file(outpatches, outfilepath, OUTBANDS)
    else
        write_tfrecord_file_new(outpatches, outfilepath, OUTBANDS)
    end

end

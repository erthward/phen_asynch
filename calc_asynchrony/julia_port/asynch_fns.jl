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

using NearestNeighbors
using Distances
using TFRecord
using Images
using Colors
using PyPlot
using Plots
using Dates
using JSON
using Glob
using GLM


#-----------
# set params
#-----------

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
if splitpath(pwd())[2] == "home"
    const ABS_DATA_DIR = "/home/drew/Desktop/stuff/berk/research/projects/seasonality/GEE_output"
else
    const ABS_DATA_DIR = "/global/home/users/drewhart/seasonality/GEE_output/SIF/"
end

"""
directory where the TFRecord data and mixerfile live
"""
const DATA_DIR = relpath(ABS_DATA_DIR)

"""
pattern that occurs just before the file number
"""
const PATT_B4_FILENUM = "Amer-"
#const PATT_B4_FILENUM = "SIF-"

"""
kernel size used by GEE to output the TFRecord files
"""
const KERNEL_SIZE = 60

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
"""
const INBANDS = ["constant", "sin_1", "cos_1", "sin_2", "cos_2"]

"""
names of the bands to save in the output TFRecord files
"""
const OUTBANDS = ["asynch", "asynch_R2", "asynch_euc", "asynch_euc_R2", "asynch_n"]

"""
max distance out to which to find and include neighbors in
each pixel's asynchrony calculation (in meters)
"""
const NEIGH_RAD = 150_000

"""
minimum distance represented by 1 degree of longitude
i.e. the distance of 1 deg lat at my highest latitude, which is 60 N/S)
"""
const MIN_DIST_°_LON=111_320

"""
a multiplicative factor to use on the first pass of filtering for candidate
neighbor pixels (before feeding those candidates into the BallTree query);
multiplying by this before subsetting using a Euclidean distance calculation
on °lat-lon values will ensure that I cast a net wide enough to include all
pixels within the nhood radius in meters
"""
const °_NEIGH_RAD_BUFF_FACTOR = 1.05

"""
use the above vars to determine thmax distance in degrees
(for the first-pass, Euclidean distance-based neighbor filtering)
"""
const MAX_DIST_° = (NEIGH_RAD/MIN_DIST_°_LON) * °_NEIGH_RAD_BUFF_FACTOR

"""
approximate radius of the earth (for the Haversine distance-based,
second-pass neighbor filtering)
NOTE: 6371008 m is a good approximation of any of the three common methods for
      representing the Earth by its 'mean' radius
"""
const EARTH_RAD = 6_371_008

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


#-----------------
# define functions
#-----------------
#
"""
Uses the data-directory path provided to get and return
lists of the input and output files' paths.
"""
function get_infile_outfile_paths(data_dir)
    # set up IO paths
    infilepaths = glob("*tfrecord", data_dir) 
    # order the infilepaths
    sort!(infilepaths)
    # exclude any previously generated output files
    infilepaths = [fp for fp in infilepaths if !occursin("_OUT", fp)]
    # get the corresponding outfilepaths
    outfilepaths = [replace(fp, r"(?<=\d{5})\.(?=tfrecord)" => "_OUT.") for fp in infilepaths]
    return infilepaths, outfilepaths
end


"""
Gets the info out of the mixerfile located in the data_dir.
"""
function read_mixer_file(data_dir)
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
function get_mixer_info(mixer_content)
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
    xmax = xmin + (dims[1] * xres)
    ymax = ymin + (dims[2] * yres)
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
function calc_patch_dimensions(mixer_content, kernel_size)
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
into a TFRecordDataset. Returns a generator that yields each next example
i.e. patch) in that dataset, parsed into an n_bands x lon x lat numpy array.
"""
function read_tfrecord_file(infilepath::String, dims, bands; default_val=DEFAULT_VAL)
    # empty Dict to hold the arrays
    arrays = Dict()
    # read the TFRecord file
    reader = TFRecordReader(infilepath)
    # beginning of example printout:
    # Example(Features(Dict{AbstractString,Feature}("constant" => Feature(#undef, FloatList(Float32[
    # loop over the examples
    for (ex_n, example) in enumerate(reader)
        # create an empty Array to hold the data
        arr = Array{Float32}(undef, 5, dims[1], dims[2])
        for (i, band) in enumerate(bands)
            band_vals = example.features.feature[band].float_list.value[1:prod(dims)]
            band_arr = reshape(band_vals, dims) |> transpose
            # replace the default missing-data val
            replace!(band_arr, default_val=>NaN)
            arr[i,:,:] = band_arr
        end
        arrays[ex_n] = arr 
    end
    # sort the dict by key (to avoid write files out with patches in diff order)
    return sort(arrays)
end


"""
Writes the output patches (i.e. rasters) to disk as a TFRecord file,
using the given filepath.
"""
function write_tfrecord_file(patches, filepath, bands)
    # instantiate the TFRecord writer
    writer = TFRecordWriter(filepath)
    for (patch_i, patch) in patches
        # create a Dict of the outbands and their arrays
        outdict = Dict()
        for (band_i, band) in enumerate(bands)
            # transpose, set all missing back to the missing-data default val,
            # then recast as Float32 and cast as a vector
            outpatch = vec(Float32.(replace(transpose(patch[band_i, :, :]),
                                            NaN=>DEFAULT_VAL)))
            outdict[band] = outpatch
        end
        # serialize to a TF example
        write(writer, outdict)
    end
    close(writer)
end


"""
Takes the overall x and y min values of the TFRecord dataset,
the x and y resolutions, the patch dimensions, and the column and row
indices of the current patch. Returns a meshgrid of the cell centers of the
current patch.
"""
function get_patch_lons_lats(xmin, ymin, xres, yres, dims, patch_j, patch_i)
    # calculate the x and y mins of the current patch
    patch_xmin = xmin + (patch_j * dims[1] * xres)
    patch_ymin = ymin + (patch_i * dims[2] * yres)

    # get lists of xs and ys of all the current patch's pixels
    xs = LinRange(patch_xmin,
                  # NOTE: start at center of leftmost pixel,
                  # step to middle of rightmost pixel
                  # (i.e. step dims[i]-1 pixels to the right),
                  # getting a list of pixel-center
                  # coordinates of total length dims[i]
                  patch_xmin + (xres * (dims[1] - 1)),
                  dims[1])
    ys = LinRange(patch_ymin,
                  patch_ymin + (yres * (dims[2] - 1)),
                  dims[2])

    # get the meshgrid of those coordinates
    gridx = xs' .* ones(length(xs))
    gridy = ys .* ones(length(ys))'
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
function make_design_matrix()
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
function calc_time_series(patch, i, j, design_mat)
    # multiply the pixel's set of coefficients by the design mat, then sum
    # all the regression terms
    # NOTE: the coeffs are array of shape (5,);
    #       the design matrix is array of shape (365, 5);
    #       the multiplication operates row by row, and elementwise on each
    #       row, returning a new (365, 5) array that, when summed across the
    #       columns (i.e. the regression's linear terms) gives a (365,) vector
    #       of fitted daily values for the pixel
    ts = vec(sum(patch[:, i, j]' .* design_mat, dims=2))
    return ts
end


"""
Calculates the simple linear regression of y on x,
where variables are column names in DataFrame df.

Returns the intercept, slope, and R^2 in a Dict.
"""
function run_linear_regression(x, y; fit_intercept=true)
    # make the design matrix
    if fit_intercept
        X = hcat(ones(length(x)), x)
    else
        X = reshape(x, length(x), 1)
    end

    # fit a model
    mod = fit(LinearModel, X, y)

    # get the R2 and coeff values
    R2 = r2(mod)
    if fit_intercept
        int, slope = coef(mod)
    else
        int = NaN
        slope = coef(mod)[1]
    end
    return Dict("slope" => slope,
                "intercept" => int, 
                "R2" => R2)
end


"""
Standardizes a 1d Array.

Returns the standardized array of identical shape.
"""
standardize(a1) = (a1.-mean(a1))/std(a1)


"""
Finds all pixel-center coordinates within neigh_dist (in meters) of
the input lat,long coordinates.

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
function get_neighbors_info(i, j, vec_is, vec_js,
                            foc_y, foc_x, vec_ys, vec_xs,
                            yres, xres, patch_dims, hav_fn, tree;
                            neigh_rad=NEIGH_RAD)
    # get the tree's data-row numbers for all neighbors
    # within the required distance
    neighs = inrange(tree, [foc_x, foc_y], NEIGH_RAD)
    #neighs = inrange(tree, [foc_x, foc_y], MAX_DIST_°)
    # use the row numbers to subset the pixels' centerpoint coords and
    # their index coordinate pairs
    neigh_xs = vec_xs[neighs]
    neigh_ys = vec_ys[neighs]
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
    #       present an issue if/when I'm crunching the ~500 m-res NIRv data...)
    #
    neigh_dists = hav_fn.(zip(neigh_xs, neigh_ys))

    # combine into a dict
    out = Dict(zip(zip(neigh_xs, neigh_ys),
                   zip(neigh_dists, zip(neigh_is, neigh_js))))
    # TODO: DELETE ME: then filter for only pixels within the nhood radius
    # TODO: DELTE ME: filter!(kv -> kv[2][1] < NEIGH_RAD, out)
    return out
end


"""
Makes vectors of the i and j indices of the patch's pixels.
The structure and order of these vectors makes them
subsettable by the output from a neighbor search on a BallTree.
"""
function get_is_js(dims)
    is = vec(repeat([1:dims[1];], dims[2]))
    js = vec(repeat([1:dims[2];], inner=dims[1]))
    return is, js
end


"""
Read the file's data from a TFRecord file, as a set of Arrays,
one per patch.

Return both that read data and a set of identically structured
output Arrays.
"""
function get_inpatches_outpatches(infilepath, inbands, dims)
    # read the file's data in as a set of examples (i.e. patches)
    inpatches = read_tfrecord_file(infilepath, dims, inbands)
    # create a data container for all the asynch patches
    outpatches = Dict()
    for (i, patch) in inpatches
        # create the output arrays, to be filled in (start as all NaNs)
        outpatches[i] = similar(patch) * NaN
    end
    return (inpatches, sort(outpatches))
end


"""
Calculates and returns the asynchrony value for a single pixel
located at position i,j in the patch.
"""
function calc_asynch_one_pixel(i, j, vec_is, vec_js,
                               foc_y, foc_x, vec_ys, vec_xs,
                               patch, outpatch, patch_n,
                               yres, xres, dims, design_mat, tree;
                               timeit=true, verbose=true)
    if verbose && timeit
        # get start time
        start = now()
    end

    if verbose
        println("\nPROCESSING PATCH $patch_n: PIXEL ($i, $j)...")
        println("\tLAT: $foc_y; LON: $foc_x")
    end

    # create lists of R2 and dist values
    R2s::Array{Float64,1} = []
    ts_dists::Array{Float64,1} = []
    geo_dists::Array{Float64,1} = []

    # calculate the focal pixel's time series (and its standardized time series)
    ts_foc = calc_time_series(patch, i, j, design_mat)
    stand_ts_foc = standardize(ts_foc)

    # make the Haversine instance and function for this cell
    # TODO: FIX THIS LINE
    hav_inst = Haversine(EARTH_RAD)
    hav_fn = function(coord_pair)
        return hav_inst((foc_x, foc_y), coord_pair)
    end

    # get the coords, dists, and array-indices
    # of all of the focal pixel's neighbors
    coords_dists_inds = get_neighbors_info(i, j, vec_is, vec_js,
                                           foc_y, foc_x, vec_ys, vec_xs,
                                           yres, xres, dims, hav_fn,
                                           tree)

    # loop over neighbors
    for (neigh_coords, neigh_info) in coords_dists_inds
        # unpack the neighbor info
        neigh_dist, neigh_inds = neigh_info
        ni, nj = neigh_inds

        # get the neighbor's time series
        ts_neigh = calc_time_series(patch, ni, nj, design_mat)
        stand_ts_neigh = standardize(ts_neigh)

        # drop this pixel if it returns NAs
        if any(isnan.(ts_neigh))
           # do nothing 
        else
            # append the distance
            append!(geo_dists, neigh_dist)

            # calculate and append the R2
            R2 = run_linear_regression(ts_foc, ts_neigh)["R2"]
            append!(R2s, R2)

            # calculate and append the Euclidean distance
            # between standardized time series
            ts_dist = euclidean(stand_ts_foc, stand_ts_neigh)
            append!(ts_dists, ts_dist)
        end
    end

    # get the slope of the overall regression of R2s on geo dist
    # NOTE: setting fit_intercept to false and subtracting 1
    #       from the array of R2s effectively
    #       fixes the intercept at R2=1
    res = run_linear_regression(geo_dists, R2s .- 1, fit_intercept=false)
    # get the slope of the overall regression of Euclidean ts dists on geo dist
    # NOTE: just setting fit_intercept to false fits ts_dist to 0 at geo_dist=0
    res_euc = run_linear_regression(geo_dists, ts_dists, fit_intercept=false)
    # extract both results into vars
    asynch = abs(res["slope"])
    asynch_R2 = res["R2"]
    asynch_euc = res_euc["slope"]
    asynch_euc_R2 = res_euc["R2"]
    # extract sample size (i.e. number of neighbors) for this focal pixel,
    # checking that there was no issue with uneven numbers of results
    @assert(length(geo_dists) == length(ts_dists) == length(R2s), ("Lengths of " *
                                                                   "geo_dists, " *
                                                                   "ts_dists, " *
                                                                   "and R2s are " *
                                                                   "not equal!"))
    asynch_n = length(geo_dists)

    # update the output patch with the new values
    outpatch[:, i, j] = [asynch, asynch_R2, asynch_euc, asynch_euc_R2, asynch_n]

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
function calc_asynch(inpatches, outpatches,
                     patch_is, patch_js, patch_ns, vec_is, vec_js,
                     xmin, ymin, xres, yres, dims, design_mat, kernel_size;
                     trim_margin=false, verbose=true, timeit=true)
    # loop over the patches from the current file
    # NOTE: each patch is of shape (n_bands, lat, lon)
    for (patch_idx, inpatch) in inpatches
        # grab important variables
        outpatch = outpatches[patch_idx]
        patch_i = patch_is[patch_idx]
        patch_j = patch_js[patch_idx]
        patch_n = patch_ns[patch_idx]

        # calculate half kernel width
        hkw = floor(Int, kernel_size/2)

        # get the lons and lats of the pixels in the current example's patch
        # as well as an array containing columns for each of the coord pairs of the pixels
        xs, ys = get_patch_lons_lats(xmin, ymin, xres, yres, dims, patch_j, patch_i)
        vec_ys = vec(ys)
        vec_xs = vec(xs)

        # get the BallTree for these lons and lats
        # make a KDTree out of all the coordinates in a patch,
        # to be used in neighbor-searching
        # TODO: DELETE ME: tree = KDTree([vec_xs vec_ys]')
        # TODO: DELETE ME: tree = BallTree([vec_xs vec_ys]')
        tree = BallTree([vec_xs vec_ys]', Haversine(EARTH_RAD))

        #----------------
        # run calculation
        #----------------

        # loop over pixels (excluding those in the kernel's margin,
        # since there's no use wasting time calculating for those)
        for i in (hkw + 1):(size(inpatch)[2]-hkw)
            for j in (hkw + 1):(size(inpatch)[3]-hkw)

                # leave the asynch output val as a NaN if the coeffs for this
                # pixel contain NaNs
                if any(isnan.(inpatch[:, i, j]))
                    if verbose
                        println("\nPROCESSING PATCH $patch_n: PIXEL ($i, $j)...")
                        lat = vec_ys[i]
                        lon = vec_xs[j]
                        println("\tLAT: $lat; LON: $lon")
                        println("\n\truntime: ---")
                        println("==================================")
                    end
                else
                    foc_y, foc_x = (ys[i,j], xs[i,j])
                    calc_asynch_one_pixel(i, j, vec_is, vec_js,
                                          foc_y, foc_x, vec_ys, vec_xs,
                                          inpatch, outpatch, patch_n,
                                          yres, xres, dims, design_mat, tree,
                                          verbose=verbose, timeit=timeit)
                end
            end
        end
    end

    # trim the half-kernel-width margin, if requested
    if trim_margin
        if verbose
            println("\n\nTrimming $hkw-cell margin around each patch.\n\n")
        end
        # subset
        trimmed_outpatches = [op[:, hkw+1:-hkw, hkw+1:-hkw] for op in outpatches]
        return trimmed_outpatches
    else
        return outpatches
    end
end


"""
Return an output Dict containing the row, column, and patch numbers,
and outfile paths (as values, each of which are sub-Dicts)
for all files (keys).
"""
function get_row_col_patch_ns_allfiles(data_dir, patt_b4_filenum)
    # set the starting row, column, and patch counters
    patch_i = 0
    patch_j = 0
    patch_n = 0

    # get the mixer file info
    mix = read_mixer_file(DATA_DIR)
    (dims, crs, xmin, ymin, xres, yres,
     patches_per_row, tot_patches) = get_mixer_info(mix)

    # get all the input and output file paths
    infilepaths, outfilepaths = get_infile_outfile_paths(DATA_DIR)

    # get the regex pattern
    patt = Regex("(?<=$PATT_B4_FILENUM)\\d{5}")

    # assert that both lists are sorted in ascending numerical order
    # NOTE: if this is not true then my method for tracking the row, col, and
    # patch numbers will break!
    for filepaths in [infilepaths, outfilepaths]
        filenums = map(x -> parse(Int32, match(patt, splitdir(x)[2]).match),
                       filepaths)
        filenums_plus1 = [1:length(filepaths);]
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
        file_dict["patch_is"] = []
        file_dict["patch_js"] = []
        file_dict["patch_ns"] = []

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


#----------------------
# get the design matrix
#----------------------

"""
Design matrix from the harmonic regression,
to be paired with the coefficients in order to
calculate the 365-day fitted time series
"""
const DESIGN_MAT = make_design_matrix()

#---------------
# get mixer info
#---------------

"""
Information read in from the JSON mixer file.
"""
mix = read_mixer_file(DATA_DIR)
const (DIMS, CRS, XMIN, YMIN, XRES, YRES,
       PATCHES_PER_ROW, TOT_PATCHES) = get_mixer_info(mix)

"""
vectors of pixel indices in the i and j array dims
(which will be subsetted by BallTree neighbor-query output)
"""
const VEC_IS, VEC_JS = get_is_js(DIMS)

#------------------------------------
# get files' row, col, and patch info
#------------------------------------

"""
Dict containing input filenames as keys and the files' patch, column, and row
numbers as values.
Needed in order to parallelize the computation across files while still
calculating lats and lons correctly for each file's pixels
"""
const FILES_DICT = get_row_col_patch_ns_allfiles(DATA_DIR, PATT_B4_FILENUM)

"""
Array containing the input filenames
"""
const FILENAMES = [fn for fn in FILES_DICT]


#= 
#----------------------------------------------------------
# define the main function, to be mapped over a worker pool
#----------------------------------------------------------
def main_fn(file_info, verbose=VERBOSE, trim_margin=TRIM_MARGIN):
    """
    Main function to be mapped over a Pool instance

    Takes the file_info for an input file, calculates
    asynchrony for that file, writes it to file, returns None.

    The argument 'file_info' must be a dict item object of the form:
      (infilepath,
        {'outfilepath': outfilepath,
         'patch_is' = [patch_i_1, patch_i_2, ..., patch_i_I],
         'patch_js' = [patch_j_1, patch_j_2, ..., patch_j_J],
         'patch_ns' = [patch_n_1, patch_n_2, ..., patch_n_N]
         }
      )
    """

    infilepath = [*file_info][0]
    print(('\nstarting job for file "%s" '
           'in process number %s') % (os.path.split(infilepath)[1],
                                          str(os.getpid())))
    file_dict = [*file_info][1]
    outfilepath = file_dict['outfilepath']
    patch_is = file_dict['patch_is']
    patch_js = file_dict['patch_js']
    patch_ns = file_dict['patch_ns']

    # read the data in, and set up the output patches' data structure
    # NOTE: by setting this up outside the calc_asynch function, 
    #       I make it so that I can run the script for a short amount of
    #       time, then retain the partial result
    #       for interactive introspection
    inpatches, outpatches = get_inpatches_outpatches(infilepath, INBANDS,
                                                     DIMS)

    print(('RUNNING ASYNCH CALC FOR FILE: '
           '"%s"') % os.path.split(infilepath)[1])
    for patch_i, patch_j, patch_n in zip(patch_is, patch_js, patch_ns):
        print('\tPATCH: %i (patch row %i, patch col %i)' % (patch_n, patch_i, patch_j))


    # run the asynchrony calculation
    outpatches = calc_asynch(inpatches, outpatches, patch_is, patch_js, patch_ns,
                             DIMS, XMIN, YMIN, XRES, YRES, DESIGN_MAT, INDICES,
                             KERNEL_SIZE, trim_margin, verbose, TIMEIT)
    print([p.shape for p in outpatches])

    # write out the asynch data
    write_tfrecord_file(outpatches, outfilepath, OUTBANDS)

=#

"""
Make a GeoArray (multilayer raster object) of the patch provided.
Pull the necessary georeferencing info from dims and the mixer_info provided.
NOTE: dims fed in separately because dims indicated in GEE's mixer file are incorrect.
"""
function make_geoarray(patch, dims, mixer_info)
    # make a GeoArray
    # NOTE: not sure why I need to transpose the patch...
    ga = GeoArray(transpose(patch))
    # get the projection info
    proj = mixer_info["projection"]
    # get and set the CRS 
    crs = parse(Int, split(proj["crs"], ':')[2]) 
    epsg!(ga, crs)
    # get and set the bounding box
    aff = proj["affine"]["doubleMatrix"]
    min_x = aff[3]
    min_y = aff[6]
    max_x = min_x + aff[1]*dims[1]
    max_y = min_y + aff[5]*dims[2]
    bbox!(ga, (min_x=min_x,
               min_y=min_y,
               max_x=max_x,
               max_y=max_y))
    return ga
end


"""
Plot the asynch, its R2s, and its sample sizes.
"""
function plot_results(outpatch)
    # plot asynch
    p1 = plot(make_geoarray(outpatch[1]), title="asynch (corr)")
    p2 = plot(make_geoarray(outpatch[2]), title="asynch (corr): R2s")
    p3 = plot(make_geoarray(outpatch[3]), title="asynch (euc)")
    p4 = plot(make_geoarray(outpatch[4]), title="asynch (euc): R2s")
    p5 = plot(make_geoarray(outpatch[5]), title="n")
    plot(p1, p2, p3, p4, p5, layout=(3,2))
end


###########################################

fp = "../../GEE_output/TEST_coeff_output_SIF_C_S_Amer-00000.tfrecord"
dims = (316, 316)
out = read_tfrecord_file(fp, dims, INBANDS)



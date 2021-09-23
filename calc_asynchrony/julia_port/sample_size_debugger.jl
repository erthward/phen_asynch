using OrderedCollections
using NearestNeighbors
using StaticArrays
using Distributed
using Statistics
using Distances
using GeoArrays
using Printf
using Images
using Colors
using Plots
using Dates
using Glob


function run_it(fn::String)::OrderedDict{Int64, Array{Float32, 3}}
    fi = (fn, FILES_DICT[fn])
    main_fn(fi)
    outdims = (316-60, 316-60)
    res, bad = get_inpatches_outpatches(replace(fn, r"-(?=\d{5})" => "-OUT-"),
                                        OUTBANDS, outdims)
    return res
end

function plot_it(res)
    p1 = plot(GeoArray(res[1][:,:,5]))
    p2 = plot(GeoArray(res[2][:,:,5]))
    plot(p1, p2)
end



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
function get_neighbors_info(i, j, 
                            vec_is, vec_js,
                            foc_y, foc_x,
                            vec_ys, vec_xs,
                            yres, xres,
                            patch_dims, hav_fn,
                            tree, neigh_rad)::Dict{Tuple{Float64, Float64}, Tuple{Float64, Tuple{Int64, Int64}}}
    # get the tree's data-row numbers for all neighbors
    # within the required distance
    neighs = inrange(tree, [foc_x, foc_y], NEIGH_RAD)
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
    return out
end





########################################




NEIGH_RAD = 150_000

EARTH_RAD = 6_371_008

DIM = 164


vec_is = vec(repeat([1:DIM;], DIM))
vec_js = vec(repeat([1:DIM;], inner=DIM))

patch_ymin = -40
patch_ymax = 40
patch_xmin = -119
patch_xmax = -79
ys = LinRange(patch_ymin, patch_ymax, DIM)
ys = ys .* ones(length(ys))'
xs = LinRange(patch_xmin, patch_xmax, DIM)
xs = xs' .* ones(length(xs))
vec_ys = vec(ys)
vec_xs = vec(xs)

yres=(patch_ymax-patch_ymin)/DIM
xres=(patch_xmax-patch_xmin)/DIM

dims = [DIM DIM]

hav_inst = Haversine(EARTH_RAD)

tree = BallTree([vec_xs vec_ys]', hav_inst)



results = zeros(DIM, DIM)
for ci in CartesianIndices(results)
    i,j = Tuple(ci)
    foc_y = ys[i, j]
    foc_x = xs[i, j]
    #println(foc_y, "   ", foc_x)
    #println(i, "   ", j)

    hav_fn = function(coord_pair)
            return hav_inst((foc_x, foc_y), coord_pair)
    end

    neighs = get_neighbors_info(i, j, vec_is, vec_js, foc_y, foc_x, vec_ys, vec_xs, yres, xres, dims, hav_fn, tree, NEIGH_RAD)
    results[i, j] = length(neighs)
    #println(length(neighs))
    #println("\n==================\n")
end

pyplot()
heatmap(results)

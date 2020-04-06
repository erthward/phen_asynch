# TODO:
    # 1. decide language

    # 2. decide if topographic correction needed

    # 3. may need to change approach for reading and extracting data, because
    # we may not be able to/want to read into memory whole time series (e.g. SIF
    # data) and/or whole global extent (e.g. SRTM)

    # 4. decide how to express and determine max neighbor-distance for
    # asynchrony analysis

    # 5. decide which topographic 'complexity' metric to use


#############################################################
########  IO                                         ########
#############################################################


# functions to read various filetypes
def read_netCDF(filename, var_names):
    # call function(s) from necessary package(s) to read data

    # subset the necessary variables

    return

def read_GeoTiff(filename):
    # call function(s) from necessary package(s) to read data
    return


# other read functions?
# .
# .
# .




# functions to read particular datasets
def read_OCO2_gridded(filename):
    data = read_netCDF(filename, var_names)
    return data


def read_OCO2_original(filename):
    data = read_netCDF(filename, var_names)
    return data


def read_TROPOMI(filename):
    data = read_netCDF(filename, var_names)
    return data


def read_ESA_CCI_LC(filename):
    data = read_GeoTiff(filename)
    return data


def read_SRTM(filename):
    data = read_GeoTiff(filename)
    return data


# NOTE: Don't know what files and formats we'll have for precip, 
# temperature, and cloud-cover yet
def read_tmp(filename):
    return

def read_ppt(filename):
    return

def read_cld(filename):
    return




#############################################################
########  SPATIAL                                    ########
#############################################################

def reproject_raster(rast_stack):
    # call function(s) from necessary package(s)
    return


def calc_topo_complexity(DEM):
    # after the metric is decided, either find a package to
    # calculate it or write up an algorithm to calculate it
    return


def get_orbital_gap_points(OCO2_orig):
    # NOTE: decide how to 'look through' all the OCO2 original data and return
    # a list of cells without coverage
    # tentative idea: create an np.array of zeros, then loop through the data
    # and flip the zeros to ones for all cells with non-null values, such that
    # the leftover zeros at the end will represent all orbital-gap cells
    return



#############################################################
########  VALIDATION                                 ########
#############################################################

def calc_gap_validation(OCO2_gridded, OCO2_orig, TROPOMI, R2_min):
    gap_points = get_orbital_gap_points(OCO2_orig)

    results = []

    for point_x, point_y in gap_points:
        # get gridded OCO2 and TROPOMI time series at the point
        OCO2_tseries = get_point_tseries(point_x, point_y, OCO2_gridded)
        TROP_tseries = get_point_tseries(point_x, point_y, TROPOMI)

        # calculate R-squared 
        R2 = calc_tseries_correlation(OCO2_tseries, TROP_tseries)

        # check threshold
        results.append(R2 >= R2_min)

    return result



#############################################################
########  VALUE EXTRACTION                           ########
#############################################################

# extract the full time series at a point (for raster stacks)
def get_point_tseries(point_x, point_y, rast_stack):
    tseries = rast_stack[point_y, point_x, :]
    return tseries

# extract the single value at a point (for raster layers)
def get_point_val(point_x, point_y, rast):
    val = rast[point_y, point_x]
    return val


#############################################################
########  NEIGHBOR EXTRACTION                        ########
#############################################################
    # TODO: devise some algorithm for finding all cells within the max
    # distance
def get_neighbors(point_x, point_y, rast, max_dist):
    # get something like a list of all neighbors,
    # expressed as tuples of (neigh_x, neigh_y, neigh_dist)
    neighs = [(), (), (), ...]
	
	# take all the cells within a certain radius of the raster point of radius size x
	rast_coord = make_array(point_x, point_y)
	rast_neighs =  st_buffer(rast_coord, dist, nQuadSegs)
    
    #filter out certain neighbors
	def sort_neighs(x)
		if x = rast_neighs_water:
            # filter out bodies of water
			return # get x and y coordinate and raster 
		else x = rast_neighs_ag:
            # filter out certain agricultural lands
			return # get x and y coordinate and raster 
   
    # Extract list of all neighbors that make the cut
	sorted_neighs = sort_neighs(x)
    
    #Extract the distance of all neighbors from ‘neigh_coord’ that make the cut
	neigh_dist = st_distance(rast_coord, sorted_neighs) 
		# Takes raster layer coordinate and neighbor coordinate 
	
	neighs = append(neigh_coord, neigh_dist)
    
	return neighs



#############################################################
######## DATA FILTERING                              ########
#############################################################

def get_tseries_lengths(rast_stack, lats, lons):
    # assuming we have something like a Numpy array that has dimensions
    # lat_steps (i) x lon_steps (j) x time_steps (k):
    lens = []
    for lat in lats:
        for lon in lons:
            lens.append(np.sum(np.invert(np.isna(rast_stack[lat, lon, :]))))

    return lens


def get_tseries_annual_coverages(rast_stack, lats, lons):
    # NOTE: need to decide how to represent and measure the proportion of the
    # calendar year that is covered by the non-null values in the time series
    # of a site
    return


def get_lc_proportions(cell_llx, cell_lly, cell_urx, cell_ury, lc):
    # NOTE: need to decide how to summarize and then filter based on the
    # proportions of different land-cover types within the footprint of a given
    # SIF grid cell (i.e. in the window bounded by (llx, lly) and (urx, ury)
    return


def get_SIF_qc_vals():
    # NOTE: don't know if this will be relevant/necessary
    return


def get_point_filter_decision(tseries_len, tseries_coverage, point_lc_props,
                              point_qc_vals):
    # make a binary decision whether or not to filter the site based on its
    # various input values
    return



#############################################################
######## FFT                                         ########
#############################################################

def calc_FFT_power(tseries, ann_freq, biann_freq):
    # use necessary pacakge to calculate the FFT
    fft_output = calc_fft(tseries)

    # discern whether this site should be considered an annual or biannual site
    # (based on which of the two frequencies has a higher power density
    ann_pow = fft_output[ann_freq]
    biann_pow = fft_output[biann_freq]
    if ann_pow > biann_pow:
        return 1, ann_pow
    else:
        return 2, biann_pow


def calc_FFT_pval(tseries, power):
    nullseries = deepcopy.deepcopy(tseries)
    # use permutation test
    null_pows = []
    for i in range(int(1e5)):
        np.random.shuffle(nullseries)
        seasonality, null_pow = calc_FFT_power(nullseries)
        null_pows.append(null_pow)

    permut_result = mean(null_peaks >= peak)

    # use frequentist test, if we have one
    freq_result = frequentist_test(tseries, peak)

    result = permut_result and freq_result

    return result, seasonality



#############################################################
######## CORRELATION                                 ########
#############################################################

def calc_series_correlation(series_x, series_y):
    # use correlation function from some package
    coeff, R2 = corr(series_x, series_y)
    return coeff, R2


def calc_rast_correlation(rast_x, rast_y):
    # NOTE: rasters would need to have been resampled to the gridded OCO2 resolution
    # use correlation function from some package
    coeff, R2 = corr(rast_x, rast_y)
    return coeff, R2



#############################################################
######## ASYNCHRONY                                  ########
#############################################################

def calc_point_asynchrony(point_x, point_y, rast_stack, max_dist):
    # create empty data structures
    neigh_dists = []
    neigh_R2s = []
    # get the point's time series
    tseries = get_point_tseries(point_x, point_y, rast_stack)
    # get all the point's neighbors
    neighs = get_neighbors(point_x, point_y, rast_array, max_dist)
    # calculate correlation between focal site and each neighbor
    for neigh_x, neigh_y, neigh_dist in neighs:
        # save neighbor's distance
        neigh_dists.append(neigh_dist)
        # get neighbor's time series
        neigh_tseries = get_point_tseries(neigh_x, neigh_y, rast_stack)
        # save R2 between focal and neigh sites' time series
        neigh_R2s.append(calc_series_correlation(tseries, neigh_tseries)[1])

    # calculate correlation between time-series correlations' R2s and their
    # inter-neighbor distances
    asynch, R2 = calc_series_correlation(neigh_dists, neigh_R2s)
    return asynch, R2


def calc_rast_asynchrony(rast_stack, lats, lons, max_dist):
    # create an empty ixj array for the asynch results
    asynch = np.full(rast_stack.shape[:2], np.nan)
    # then loop over all cells and fill it up
    for i, lat in enumerate(lats):
        for j, lon in enumerate(lons):
            asynch[i, j] = calc_point_asynchrony(lon, lat, rast_stack,
                                                     max_dist)
    return asynch



#############################################################
######## ANALYSIS                                    ########
#############################################################

def do_analysis():
    # read in and pre-process (e.g. reconcile projections) all data

    # run orbital-gap validation on gridded data

    # filter unusable cells out of data, using:
        # 1. number of usable neighbor sites (based on the length and annual
        #        coverage of their time series)
        # 2. land cover classes
        # 3. QC flags from SIF data, if applicable
        # 4. passing of FFT significance test for annual or biannual signal 

    # calculate global raster of topographic complexity

    # calculate global SIF asynchrony map

    # calculate global asynchrony maps of precip, temp, and clouds

    # calculate raster correlations between SIF asynchrony and:
        # a. latitude and longitude
        # b. elevation and topographic complexity
        # c. climate, and asynchrony of climate

    return

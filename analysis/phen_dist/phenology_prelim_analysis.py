

### SKELETON CODE FOR PHENOLOGY DATA PRELIMINARY ANALYSIS FUNCTION ###

# IMPORTANT NOTE ABOUT gBIF DATA TO KEEP IN MIND: ALL observations have life stage = "Flowering" but this column isn;t included because it is redundant

# extract a single species from filtered gBIF dataset (goal is to do this for all species in dataset) 
    # using column name `verbatimScientificName` or column `species` or `speciesKey` which assigns a unique 7 digit code to every unique species


# extract the longitude and latitude of a single initial observation among all observations for the selected species (siteA)
    # might want to form a matrix
    # using column names `decimalLatitude` and `decimalLongitude`


# extract the date from initial observation selected above (since we know this observation was life stage = "Flowering") 
    # column name `eventDate` is formatted as so: "1974-12-27T00:00:00"
    # columns `day`, `month`, and `year` provide the numeric month day and year in the following form: "27"	"12"	"1974" 
   

# convert extracted date from initial observation into the numeric day on the calendar and then into a value on the unit circle
    def get_doy_from_month_day(year, month, day):
        return (datetime.datetime(year, month, day).timetuple().tm_yday)/365*(2*np.pi)
    
    siteA_doy_rad = get_doy_from_month_day(year, month, day)

# select another observation from the observations of a selected species (goal is to iterate over all observations for a single species) 
    # select observation whose `gbifID` is different from that of the initial observation
    # `gbifID` assigns a unique 10 digit number to all unique observations
    
  # extract the longitude and latitude of the selected observation among all observations for the selected species (siteB)
      # might want to form a matrix

  # extract the date from the selected observation
  
  # convert extracted date from selected observation into the numeric day on the calendar and then into a value on the unit circle
        siteB_doy_rad = get_doy_from_month_day(year, month, day)

# compare observation from siteA and siteB by taking the chord distance of the day of year (this is the phenological distance of the two species)
  # do this step for all pairs of observations within a species (initial observation, all other observations for that species) 
  # where "siteB_doy_rad" is a vector of all radian DOY values for all observations in a particular species not including the inital observation
    # such that the length is (number_of_species_observations - 1)
       
       chord_dist = [np.sin(np.abs(siteA_doy_rad - n)/2) for n in siteB_doy_rad]
  
        
##### NOTE: THINGS STARTED GETTING SHAKY STARTING HERE 

# extract the coefficients for the fitted time series from Earth Engine for the site locations of all observations for a single species 
# using the coordinates from above
    # get the annual frequency coefficients by substituting a vector of days into: (day_vec/365) * 360 * 2*pi
    # get the biannual frequency coefficients by substituting a vector of days into: (day_vec/(365*2)) * 360 * 2*pi


# create a 365 X 5 model matrix containing a column for day of year and the other four columns contain the sine/cosine frequencies
    # ------> create two model matrices? one for siteA and another for siteB? 


# construnct the fitted seasonal function by multiplying the coeffcients and model matrix together 
    # Do this for all pairs of observations within a species (initial observation, all other observations for that species)
    # to get set of pairwise distances for a single species


# compute the euclidean distance between the two signal functions for siteA and siteB
    # Do this for all pairs of observations within a species


# Run a regression and extract the summary statistics
    # ------> is this alinear regression? or a robust regression? 
        # DETH: this would be a Mantel test (or a multi-Mantel test, which would allow us to partial out geographic distance)
        #       this is a type of multiple-matrix regression, where the eqxn is of the normal form of a regression (Y ~ beta_1*X_1 + ... + beta_n*X_n), 
        #       but where Y, X_1 ... X_n are matrices of pairwise distances (rather than vectors of IID scalars)
        
    # ------> run one regression for the fitted seasonality for siteA and another regression for the fitted seasonality for siteB? 
        # DETH: each species would have a single regression, of the form D_phen ~ beta * D_seas + beta * D_geog, where the Ds are the pairwise distance matrices,
        # the betas are the fitted coeffs the regression estimates, and each matrix is of the form:
        #
        #        1  2  3  ...   4
        #        10
        #        2  0
        #        3   0
        #        .    
        #        .
        #        .
        #        4              0
        #
        # in which each number is a separate sample, the off-diagonals are the pairwise distances (whether phenological, seasonal, or geographic)
        # between samples, the matrix is thus symmetric, and the diag is necessarily 0
        
    # ------> is this to determine if our fitted seasonality curves fit the "observed" seasonality at each site from the gBIF data?
        # DETH: this is to test whether the map of seasonality at the vegetation stand scale, which we estimated
        #       from remotely sensed data, is a significant predcitor of the phenological sasonality of individual species observed in the field



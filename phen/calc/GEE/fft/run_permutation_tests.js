/* 
EXPORTS:
  runPermutations(harmonics, bandsForViz, imgColl, dependent, nIts, maxPVal)
    harmonics: ee.List of harmonic frequencies to include in the harmonic regression
    bandsForViz: ee.List of the regression-output bands to be used for a visualization
    imgColl: imageCollection to be filtered based on permutation testing of R^2s
    dependent: String name of the dependent variable's band
    nIts: the number of iterations to use in the permutation test
    maxPVal: max p-value (above which a pixel will be filtered out)
  

*/


//--------
// imports
//--------

var hReg = require('users/drewhart/seasonality/:fft/run_harmonic_regression.js');


//----------
// functions
//----------

// return a permuted ImageCollection
var permuteImageCollection = function(imgColl, seed){
  // get the original list of start times
  var origT = imgColl.aggregate_array('t');
  // use a random column (with n, the iteration number, as the seed) to permute the IC
  var permuted = imgColl.randomColumn('rand', seed).sort('rand');
  // get the permuted list of start times
  var permutedT = permuted.aggregate_array('t');
  // remap the permuted start times in place of the original start times
  // NOTE: this is necessary so that the actual time data that is used to
  // create the covariates is actually permuted with respect to the regressand
  var out = permuted.remap(origT, permutedT, 't');
  return out;
};


// main permutation-test function
exports.runPermutations = function(harmonics, bandsForViz, imgColl, dependent, nIts, maxPVal, seedStart){

  // run the unpermuted harmonic regression and get its R^2 Image
  var trueR2 = ee.ImageCollection(hReg.calcHarmonicRegression(imgColl, dependent,
                                                          harmonics, bandsForViz))
    .first()
    .select('R2');
  
  
  // define functions needed inside this scope:
  
  // function to permute the target ImageCollection, run harmonic regression
  // on it, and return the output R^2 Image
  var getPermR2 = function(seed){
    var permuted = permuteImageCollection(imgColl, seed);
    var out = ee.ImageCollection(hReg.calcHarmonicRegression(permuted, dependent,
                                                             harmonics, bandsForViz))
      .first()
      .select('R2');
  return out;
  };
    
  
  // Function to use for iteration over list of serial integers (of len nIts)
  var addAnotherPermutationTest = function(curr_seed, img_from_previous_seed){
    var permR2 = ee.Image(getPermR2(ee.Number(curr_seed)));
    var output = ee.Image(img_from_previous_seed).add(permR2.gte(ee.Image(trueR2)).rename('fail_counts')).cast({'fail_counts': 'uint8'});
    return output;
  };
  
    
  // Function to calculate an Image containing the
  // number of permutation-test fail counts
  // with nIts permutation iterations
  var calcFailCounts = function(nIts, seedStart){
    // Calculate a list of serial integers, one for each permutation iteration
    var seeds = ee.List.sequence(seedStart, ee.Number(seedStart).add(ee.Number(nIts)).subtract(1));
    // create the starting Image of fail counts (started at 0, as no counts have yet accumulated)
    var firstImage = ee.Image.constant(0)
      // rename the band to match the band name used in each true-permuted R2 comparison
      .rename('fail_counts')
      // cast the band to an unsigned 8-bit integer
      // (minimizes memory usage, just need to be sure we don't exceeed 232 total permutation tests)
      .cast({'fail_counts': 'uint8'});
    // get the fail-count image
    var failCounts = ee.Image(seeds.iterate(addAnotherPermutationTest, firstImage))
      .rename('fail_counts');
    return failCounts;
  };
  
  
  // calculate the fail count Image
  var failCounts = calcFailCounts(nIts, seedStart)
  // reproject, to make sure that the whole series of calculations happens in the input data's CRS
    .reproject(imgColl.first().projection());
  return failCounts;  
};
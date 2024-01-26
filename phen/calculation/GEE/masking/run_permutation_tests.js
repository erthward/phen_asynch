/*

TODO:

- any options for relieving memory limitations??
  (given that I use no reduce or sample functions
  that have parallelScale and tileScale args)
  
*/


/* 
EXPORTS:
  runPermutationTest(imgColl, dependent, nIts, maxPVal)
    imgColl: imageCollection to be filtered based on permutation testing of R^2s
    dependent: String name of the dependent variable's band
    nIts: the number of iterations to use in the permutation test
    maxPVal: max p-value (above which a pixel will be filtered out)
  

*/



//////////////////////////////
// LOAD AND PREP DATASETS, ETC
//////////////////////////////


//-------
// params
//-------

var nIts = 1000;

var export_it = false;

//--------
// imports
//--------

var hReg = require('users/drewhart/seasonality/:fft/harmonic_regression.js');


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
exports.runPermutations = function(imgColl, dependent, nIts, maxPVal, seedStart){
  
  // set the harmonic regression arguments that don't need to be varied when runPermutations is called
  var harmonics = ee.List([1, 2]);
  var bandsForViz = ee.List(['sin_1', 'cos_1']);
  
  // now run the unpermuted harmonic regression and get its R^2 Image
  var R2 = ee.ImageCollection(hReg.calcHarmonicRegression(imgColl, dependent,
                                                          harmonics, bandsForViz))
    .first()
    .select('R2');
  
  
  // define functions needed inside this scope:
  
  // function to permute the target ImageCollection, run harmonic regression
  // on it, and return the output R^2 Image
  var getHarmRegR2 = function(seed){
    var permuted = permuteImageCollection(imgColl, seed);
    var out = ee.ImageCollection(hReg.calcHarmonicRegression(permuted, dependent,
                                                             harmonics, bandsForViz))
      .first()
      .select('R2');
  return out;
  };
    
  
  // Function to use for iteration over list of serial integers (of len nIts)
  var addAnotherPermutationTest = function(curr_seed, previous){
    var permutTestRes = ee.Image(getHarmRegR2(ee.Number(curr_seed)));
    var output = ee.Image(previous).add(permutTestRes.gte(ee.Image(R2)));
    return output;
  };
  
    
  // Function to calculate significance of an image's values
  // (for the only band in that image)
  // using a permutation-based significance test
  // with nIts permutation iterations
  var calcFailCounts = function(img, nIts, seedStart){
    // get the band name
    //var band = ee.String(ee.List(img.bandNames()).get(0));
    // Calculate a list of serial integers, one for each permutation iteration
    var ntimes = ee.List.sequence(seedStart, ee.Number(seedStart).add(ee.Number(nIts)).subtract(1));
    // create the starting Image
    var firstImage = ee.Image.constant(0)
      // rename the band
      .rename('R2')
      // cast the band to a 64-but integer
      .cast({'R2': 'long'});
    // get the p-value image
    var failCounts = ee.Image(ntimes.iterate(addAnotherPermutationTest, firstImage))
      //.divide(nIts)
      //.rename('p_value');
      .rename('fail_counts');
    //return p;
    return failCounts;
  };
  
  
  // then calculate the empirical p-vals from that collection of R2s,
  // and return mask
  var failCounts = calcFailCounts(R2, nIts, seedStart)
  // mask pixels with p-val > threshold,
  //  .lte(maxPVal)
  // reproject, to make sure that the whole series of calculations happens in the input data's CRS
    .reproject(imgColl.first().projection());
  return failCounts;  
};
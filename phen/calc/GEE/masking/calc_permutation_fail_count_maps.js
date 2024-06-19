// calculates a map of the number of times each cell's R2 in a randomly permuted SIF time series
// is greater than that cell's R2 in the unpermuted (i.e., true) time series
// (i.e., the number of times a cell 'fails' the permutation test)
// results will be added together then divided by the total
// number of tests run to develop a global p-value map, for QC data masking


//var metric = 'SIF';
var metric = 'NIRv';


// NOTE: DETH: 07/01/22: FOR SIF, I RAN THIS SCRIPT 4 TIMES, EACH TIME FOR A BATCH OF 25 TESTS
//                       AND USING 0, 25, 50, OR 75 AS THE seedStart VALUES, IN ORDER TO PRODUCE
//                       THE 4 MAPS THAT I THEN COMBINED INTO THE PERMUTATION-TEST SIGNIFICANCE MAP

// NOTE: DETH: 08/03/22: FOR NIRv, I RAN THIS SCRIPT ONCE EACH FOR SEEDS 0-4, 5-9, 10-14, 15-20, & 20-25,
//                      THEN COMBINED INTO THE PERMUTATION-TEST SIGNIFICANCE MAP




// number of harmonics to use in harmonic regression (2 means one annual, one semiannual)
var nFreq = 2;
var harmonics = ee.List.sequence(1, nFreq); 
// and bands to use for calculating phase and amplitude
var bandsForViz = ee.List(['sin_1', 'cos_1']);

// number of permutations,
// max p-value beyond which to mask out a pixel (i.e., alpha value)
// and number to use to start seq of integers that serve as seeds for each iteration
var nIts = 5;
var maxPVal = 0.01;
//var seedStart = 0;
//var seedStart = 25;
//var seedStart = 50;
var seedStart = 0;

// requisite functions for...
// data-loading
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
// significance-testing
var runPerms = require('users/drewhart/seasonality/:masking/run_permutation_tests.js');
// harmonic regression
var hReg = require('users/drewhart/seasonality/:fft/harmonic_regression.js');


// read in the SIF dataset and get its band name, projection, and scale
var minWaterOccur = 100; // NOTE: don't mask out water for significance tests
if (metric == 'SIF'){
  var dataset = gDat.getSIFData(minWaterOccur);
  var dependent = 'b1';
} else{
  var dataset = gDat.getNIRvData(minWaterOccur, 3);
  var dependent = 'NIRv';
}

var proj = dataset.first().projection();
var scale = dataset.first().projection().nominalScale();
  

// run permutation test and get significance mask mask
// NOTE: leaving permutation test function as returning a mask, not a masked ImgColl,
//       because I may wind up wanting to produce and export the mask (or multiples of
//       them) separately, then read them in, combine, them, and use them
var perms = runPerms.runPermutations(dataset, dependent, nIts, maxPVal, seedStart);

// set task and asset names
var taskAndAssetNameBase = ee.String(metric).cat('_permutations_seeds');
var taskAndAssetName = ee.String(taskAndAssetNameBase)
  .cat(ee.String(ee.Number(seedStart).toInt()))
  .cat('-')
  .cat(ee.String(ee.Number(seedStart + nIts - 1).toInt()))
  .getInfo();

// export as an asset
var roi =
  ee.Geometry.Polygon(
      [[[-165, 60],
        [-165, -60],
        [180, -60],
        [180, 60],
        [-165,60]]], null, false);
Export.image.toAsset({
  image: perms.select('fail_counts'),
  description: taskAndAssetName,
  assetId: taskAndAssetName,
  region: roi,
  scale: scale.getInfo(),
  maxPixels: 6e9,
  pyramidingPolicy: {
    '.default': 'mean',
}});
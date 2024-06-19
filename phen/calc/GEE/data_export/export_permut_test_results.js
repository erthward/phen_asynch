//-----
// TODO
//-----


//------
// NOTES
//------


// IMPORTANT CONCERNS, POTENTIAL COMPLICATIONS:

// 0. Code should be stored in appropriately named files that are saved in appropriate
//    subdirectories, so that things are easy to browse, read, and understand later on.
//    Thus, to the extent necessary, we should reorganize our files, and we can delete
//    files that are no longer useful (or move to 'scratch' if unsure but the code 
//    is not part of our main processing pipeline).

// 1. Will be super important to give careful attention to how to manage the scale
//    and projection desired for each step!!

// 2. Not all neighboring pixels will have the same time series in MODIS! Need to interpolate
// to daily? Or some other method of reconciling before running correlations?

// 3. For two time series with identical seasonality, but one with higher max or lower min
// vals than the other, the R-squared between them will be reduced in a way that does not
// actually reflect seasonal decoupling (and this will not be 'fixed' by just normalizing
// a time-series to itself). What to do about this??




//-------
// PARAMS
//-------

// max p-value to use for the significance test,
// and number to use to start seq of integers that serve as seeds for each iteration
var nIts = ee.Number(5);
var maxPVal = 0.05;
var seedStart = ee.Number(0);


// number of harmonics to use in harmonic regression ('2' means: one annual, one semiannual)
var nFreq = 2;
var harmonics = ee.List.sequence(1, nFreq); 
// and bands to use for calculating phase and amplitude
var bandsForViz = ee.List(['sin_1', 'cos_1']);


// dataset to use (String that will be later used to set the dataset ImgColl and the dependent var)
//var datasetName = 'SIF';
var datasetName = 'NIRvP';

// min percent of observ (over time) at a pix to cause it to be masked
var minWaterOccur = 80;

// min percent valid land cover among the land-cover pixels that fall inside
// each input dataset (NIRvP, SIF) pixel (pixels with less than this will be filtered out)
var minPctValidLC = 0.75;

// min perecent of non-null values needed in the time series 
// (less than this will be filtered)
var minPctTSLen = 0.5;

// number of permutations,

// max neighbor distance to use in calculating asynchrony (in meters)
var maxNeighDist = ee.Number(15000);

// min percent of valid neighbors that a pixel must have
// in order to be included in the output asynchrony map
var minPctValidNeigh = 0.5;


// output controls
var verbose = true;
var map_intermediates = false;
var map = false;
var custom_palette = true;
var export_result = true;




//---------
// REQUIRES
//---------

// data-loading
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
// filtering
var maskLC = require('users/drewhart/seasonality/:masking/mask_landcover.js');
var maskShortTS = require('users/drewhart/seasonality/:masking/mask_short_time_series.js');
var getSigMask = require('users/drewhart/seasonality/:masking/get_n_permut_test_results.js');
// harmonic regression
var hReg = require('users/drewhart/seasonality/:fft/harmonic_regression.js');
// annual fitted time series
var fitYr = require('users/drewhart/seasonality/:fft/fit_single_year_daily_values.js');
// asynchrony
var asynch = require('users/drewhart/seasonality/:asynchrony/calc_asynchrony.js');
// map stuff
var cbar = require('users/drewhart/seasonality/:viz/make_colorbar.js');
var palettes = require('users/gena/packages:palettes'); 



//--------
// IMPORTS
//--------

// read in the dataset needed
if (datasetName === 'NIRvP'){
  // NIRvP, its dependent variable band name, proj, scale, and nNeigh
  var dataset = gDat.getNIRvPData(minWaterOccur)
    // filter for a shorter date range, to try to get around memory-usage error
    // NOTE: for now taking a 10-year time slice to see how that does...
    .filterDate('2009-12-31T00:00', '2019-12-31T00:00');
  var dependent = 'NIRvP';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
  var nNeigh = maxNeighDist.divide(scale).ceil(); // num neighs for asynch calc
  
  // or SIF, its dependent variable band name, proj, scale, and nNeigh
} else if (datasetName === 'SIF'){
  var dataset = gDat.getSIFData(minWaterOccur);
  var dependent = 'b1';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
  var nNeigh = maxNeighDist.divide(scale).ceil(); // num neighs for asynch calc
}




//==================
// MASKING/FILTERING
//==================

//------------------
// filter land cover
//------------------
// NOTE: this will also reduce the temporal extent of the data
//       if the data's original temporal extent falls outside that of the LC data!
var LCMasked = maskLC.maskLandCover(dataset, minPctValidLC, proj);


//-------------------
// filter missingness
//-------------------
var LCTSMasked = maskShortTS.maskShortTimeSeries(LCMasked, dependent, minPctTSLen);


//---------------------------------------------------
// filter by significance of permutation test of R^2s
//---------------------------------------------------
// calculate regression with nFreq freqs
var reg = ee.ImageCollection(hReg.calcHarmonicRegression(LCTSMasked, dependent,
                                                         harmonics, bandsForViz));
// run permutation test and get mask
// NOTE: leaving permutation test function as returning a mask, not a masked ImgColl,
//       because I may wind up wanting to produce and export the mask (or multiples of
//       them) separately, then read them in here to use them
var testRes = getSigMask.runPermutationTest(LCTSMasked, dependent, nIts, maxPVal, seedStart);

if (export_result){
  var roi = 
      ee.Geometry.Polygon(
          [[[-179, 60],
            [-179, -10],
            [-34, -10],
            [-34, 60],
            [-179,60]]], null, false);
  
  // Export the image to an Earth Engine asset.
  var taskAndAssetName = ee.String('testRes_')
    .cat(datasetName)
    .cat(ee.String('_seeds_'))
    .cat(ee.String(seedStart.toInt()))
    .cat(ee.String('_'))
    .cat(ee.String(seedStart.add(nIts).subtract(1).toInt()))
    .getInfo();
  Export.image.toAsset({
    image: testRes.select('n_failed_tests'),
    description: taskAndAssetName,
    assetId: taskAndAssetName,
    scale: scale.getInfo(),
    region: roi,
    maxPixels: 4e9,
    pyramidingPolicy: {
      '.default': 'mean',
    }
  });
}

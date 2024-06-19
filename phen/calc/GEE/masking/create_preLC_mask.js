// Make pre-LC-masking mask map, both for default and strict masking modes


//-------
// PARAMS
//-------


// . . . . . . . . . 
// masking/QA params
// . . . . . . . . .

// MASKING MODE:
var maskingMode = 'default';
//var maskingMode = 'strict';


// dataset to use (String that will be later used to set the dataset ImgColl and the dependent var)
// NOTE: JUST USING TS LENGTH OF NIRv DATA TO MASK BOTH NIRv AND SIF DATA
var datasetName = 'NIRv';



// DEFAULT MASKING PARAMS
// _______________________

// minimum Pielou evenness of MODIS time series 
var minTSEvenness = 0.8;

// minimum pct time in majority (i.e., mode) MODIS lc type
var minPctTimeLCMode = 8/10; // NOTE: 8/10 = 8 of the 10 years of LC data used must be the same LC type

// min perecent of non-null values needed in the time series 
// (less than this will be filtered)
var minPctTSLen = 0.5;

// ______________________
// min percent of observations (over time) at a pix to cause it to be masked
// NOTE: 100 IN THIS CASE, BECAUSE DO NOT WANT WATER TO FACTOR IN AT THIS STEP
var minWaterOccur = 100; // NOTE: expressed as pct of 100, not fraction of 1




//-----------------------------------
// LOAD REQUIRED FUNCTIONS AND ASSETS
//-----------------------------------

// filtering
var maskLC = require('users/drewhart/seasonality/:masking/mask_landcover.js');
var maskShortTS = require('users/drewhart/seasonality/:masking/mask_short_time_series.js');
var pctTimeLCMode = ee.Image('users/drewhart/pct_time_in_MODIS_lc_mode');
var tsEvenness = ee.Image('users/drewhart/NIRv_monthly_Pielou_evenness_10YR')
  // NOTE: some pixels are getting values erroneously >1, but I spent a long time inspecting and it appears
  //       only to happen offshore and I manually checked that the Pielou evenness code is working for
  //       real data at some pixels, so I feel confident in the result and will just clamp values to 1
  .clamp({low: 0, high: 1});

// data-loading
var gDat = require('users/drewhart/seasonality/:io/get_data.js');



//--------
// IMPORTS
//--------

// read in the dataset needed
if (datasetName === 'NIRv'){
  // NIRv, its dependent variable band name, proj, scale
  var dataset = gDat.getNIRvData(minWaterOccur, 10);
  var dependent = 'NIRv';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
    
  // or SIF, its dependent variable band name, proj, scale
} else if (datasetName === 'SIF'){
  var dataset = gDat.getSIFData(minWaterOccur);
  var dependent = 'b1';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
}



//===========
// MASKING/QA
//===========

//--------------------------------------------
// mask for even time series and for stable LC
//--------------------------------------------
var evennessMask = tsEvenness.gte(minTSEvenness);
var lcStabilityMask = pctTimeLCMode.gte(minPctTimeLCMode);
  
  
//-------------------
// filter missingness
//-------------------
var shortTSMask = maskShortTS.maskShortTimeSeries(dataset, dependent, minPctTSLen);  
 
 
//---------------------------------------------------
// filter by significance of permutation test of R^2s
//---------------------------------------------------
// NOTE: USE NIRv PERMUTATION TESTS FOR DEFAULT MASKING, OR SIF FOR STRICT
if (maskingMode == 'default'){
  var signifMask = ee.Image('users/drewhart/NIRv_permutation_significance_mask');
} else if (maskingMode == 'strict'){
  var signifMask = ee.Image('users/drewhart/SIF_permutation_significance_mask');
}


//---------------------
// compose overall mask
//---------------------

var overallMask = evennessMask.and(lcStabilityMask).and(shortTSMask).and(signifMask);

//----------
// EXPORT IT
//----------
var roi =
  ee.Geometry.Polygon(
      [[[-165, 60],
        [-165, -60],
        [180, -60],
        [180, 60],
        [-165,60]]], null, false);
  

var scale = dataset
  .first()
  .projection()
  .nominalScale()
  .getInfo();


if (maskingMode == 'strict'){
  var taskAndAssetName = 'phen_mask_STRICT';
} else if (maskingMode == 'default'){
  var taskAndAssetName = 'phen_mask_DEFAULT';
}
  
// Export the image, specifying scale and region.

Export.image.toAsset({
image: overallMask,
description: taskAndAssetName,
assetId: taskAndAssetName,
scale: scale,
region: roi,
});

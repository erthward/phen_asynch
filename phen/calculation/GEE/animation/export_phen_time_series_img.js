// BEHAVIORAL PARAMS
////////////////////

// dataset to use
//var dataset = 'NIRv';
var dataset = 'SIF';

// cut the time series down from 365 days to something shorter?
// (just for debugging)
var cutTimeSeries = true;
// length to reduce time series to (if cutTimeSeries is true)
var tsLimit = 5;

// index number of the time series to start on when subsetting
// annual time series to ImageCollection that will be exported
var tsStart = 0;
// ensure that it is 0, if time series should not be cut
if (!cutTimeSeries){
  var tsStart = 0;
}


// IMPORT CODE
//////////////

// calculate annual fitted time series
var fitYr = require('users/drewhart/seasonality:fft/fit_single_year_daily_values.js');
// normalize time series 
var norm_fit_vals = require('users/drewhart/seasonality:fft/normalize_fitted_time_series.js');
// export image collection to asset
var ImgCollExport = require('users/drewhart/seasonality:data_export/export_imgColl_toAsset.js');



// READ IN DATA AND PREP DATA
/////////////////////////////

/// READ IN DATA
// harmonic regression coeffs for global data 
if (dataset == 'NIRv'){
  var coeffs = ee.Image('users/drewhart/NIRv_global_coeffs')
    .rename(['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']);
} else if (dataset == 'SIF'){
  var coeffs = ee.Image('users/drewhart/SIF_global_coeffs')
    .rename(['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']);
}

/// CALC FITTED TIME SERIES
/* get a one-year daily time series of fitted values from the haromic regression
  from 'calcHarmonicRegression' function for global SIF data & NIRv data
  w/ coeffs from harmonic regression */

var fittedSingleYearDVImgCollCoeffs = ee.ImageCollection(fitYr.
  fitSingleYearDailyValues(coeffs, coeffs.projection()));
 

/// NORMALIZE FITTED TIME SERIES
/* compute ImageCollection of normed fitted single year daily values from the 
  haromic regression for global SIF data & NIRv data w/ coeffs from harmonic 
  regression */
var normFitValsCoeffs = ee.ImageCollection(norm_fit_vals.
  normalizeFittedVals(fittedSingleYearDVImgCollCoeffs));
if (cutTimeSeries){
  var normFitValsCoeffs = normFitValsCoeffs.limit(tsLimit);
}

print('normFitValsCoeffs', normFitValsCoeffs);

// EXPORT GLOBAL, ANNULA TIME SERIES TO IMAGE ASSETS 
////////////////////////////////////////////////////

// export params
var globalROI = 
    ee.Geometry.Polygon(
        [[[-180, -90],
          [180, -90],
          [180, 90],
          [-180, 90],
          [-180, -90]]], null, false);
var assetId = dataset; // name of asset; will be of the form 'SIF_fitted_day_<IMG_NUM>'
var scale = normFitValsCoeffs.first().projection().nominalScale().getInfo();
var crs = normFitValsCoeffs.first().projection().crs().getInfo(); // 'EPSG:4326'
var maxPixels = 4e9;
var stepSize = 1; // return images for every n days

// collapse into a single, multi-band image
var multiBandImg = normFitValsCoeffs.toBands()
  .rename(ee.List.sequence(1, normFitValsCoeffs.size()).map(function (n) {
      return ee.String('day_').cat(ee.String(ee.Number(n).toInt()));}));
  
print('EXPORTING: ', multiBandImg);





// export each ith image as an Asset 


/////// ATTEMPT 1
// Returns this error message 7 times: "Cannot read properties of null (reading 'bza')
// Exports a file that cannot be opened 
/*
Export.image.toAsset({
  image: multiBandImg, 
  description: ee.String(assetId).cat('_ann_fit_vals').getInfo(),
  assetId: ee.String(assetId).cat('_ann_fit_vals').getInfo(),
  region: globalROI,
  scale: scale,
  crs: crs,
  maxPixels: maxPixels
  });
*/


/////// ATTEMPT 2
/* brute-force test attempt since using normFitValsCoeffs didn't work: 
  create a new image collection by piecing one together */
/*
var newImgCol = ee.ImageCollection([multiBandImg.select("day_1"), multiBandImg.select("day_2")]);

// This works
print("size", newImgCol.size().getInfo());
*/

/* BUT, the function below returns the following message:
"In users/drewhart/seasonality:data_export/export_imgColl_toAsset.js collection.size is not a function" */
/*
ImgCollExport.exportImageCollectionToAsset(
            newImgCol,
            ee.String(assetId).cat('_ann_fit_vals').getInfo(),
            globalROI,
            scale,
            crs,
            maxPixels,
            tsStart,
            stepSize
            );
*/


/////// ATTEMPT 3
/* The function below returns the following message:
"In users/drewhart/seasonality:data_export/export_imgColl_toAsset.js collection.size is not a function" */
ImgCollExport.exportImageCollectionToAsset(
            multiBandImg, // from above
            ee.String(assetId).cat('_ann_fit_vals').getInfo(),
            globalROI,
            scale,
            crs,
            maxPixels,
            tsStart,
            stepSize
            );


/* Before running the task, click to export it to the Google Drive and export as a GEO_TIFF */


                                         
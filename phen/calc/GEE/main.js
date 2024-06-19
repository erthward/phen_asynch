// This is the main script used to preprocess data for our global study of
// phenology and phenological asynchrony (Terasaki Hart et al. 2023)


// TODO:
  // 1. rerun SIF datasets with 8/10 LC stability (instead of mistaken 7/9)
  // 2. run final-version data for both phen and for all climate datasets


//-------
// PARAMS
//-------

//. . . . . . . . . . . . . . 
// harmonic regression params
//. . . . . . . . . . . . . .

// number of harmonics to use in harmonic regression (2 means one annual, one semiannual)
var nFreq = 2;
var harmonics = ee.List.sequence(1, nFreq); 
// and bands to use for calculating phase and amplitude
var bandsForViz = ee.List(['sin_1', 'cos_1']);


// dataset to use (String that will be later used to set the dataset ImgColl and the dependent var)
var datasetName = 'NIRv';
//var datasetName = 'SIF';
//var datasetName = 'MODISCloud';
//var datasetName = 'TerraClimate';

// number of years of MODIS NIRv data to use
var nYearsNIRv = 10;

// variable to use, if using TerraClimate data
//var climateVar = 'pr';
//var climateVar = 'def';
//var climateVar = 'tmmn';
var climateVar = 'tmmx';


// . . . . . . . . . 
// masking/QA params
// . . . . . . . . .

// MASKING MODE:
var maskingMode = 'default';
//var maskingMode = 'strict';


// ______________________
// min percent of observations (over time) at a pix to cause it to be masked
var minWaterOccur = 50; // NOTE: expressed as pct of 100, not fraction of 1


// STRICT MASKING PARAMS

// min percent valid land cover among the land-cover pixels that fall inside  
// each input dataset (NIRv, SIF) pixel (pixels with less than this will be filtered out)
var minPctValidLC = 0.8;




// . . . . . . . . . 
// behavioral params
// . . . . . . . . . 

// output controls
var verbose = true;
var map_intermediates = false;
var map = true;
var map_result = false;
var custom_palette = true;
var export_result = true;


// . . . . . . . 
// export params
// . . . . . . . 

var in_chunks = false;
var maxFileSize = 10000000; // in bytes
// DETH: changing from 256x256 patches to 300x300
//       because that evenly breaks up the [(-165,60), (180,-60)] quasi-global study area
//       into 184 even-sized patches (23 cols x 8 rows, 
//       with each patch being 476*476 after adding the 176-cell margins);
//       each patch should be approx. 4.3 MB, so if I set the maxFileSize to 5_000_000 bytes then
//       I should get 184 4.3MB files...
var patchSize = 300; //square-patch size, in pixels
var maxPixels = 25000000; // maxPixels argument fed to export.Image.toDrive()
var fileFormat = 'TFRecord';// 'Asset', 'GeoTIFF'
// max neighbor distance to be used in calculating asynchrony (in meters)
var maxNeighDist = ee.Number(180000);
// factor by which to multiply the calculated min-needed kernel size,
// to ensure adequate kernel size (in pixels) to really get all neighbors within max neigh dist (in m)
var kernelSizeBufferFactor = 1.8;
// abbreviated name of region enclosed by the bounding box (within which outputs will be generated)
var regionAbbrev = ee.String('global');
var folder = 'phen_outputs_from_GEE';



//-----------------------------------
// LOAD REQUIRED FUNCTIONS AND ASSETS
//-----------------------------------

// filtering
var maskLC = require('users/drewhart/seasonality/:masking/mask_landcover.js');
if (maskingMode == 'default'){
  var preLCMask = ee.Image('users/drewhart/phen_mask_DEFAULT');
} else if (maskingMode == 'strict'){
  var preLCMask = ee.Image('users/drewhart/phen_mask_STRICT');
}

// data-loading
var gDat = require('users/drewhart/seasonality/:io/get_data.js');

// harmonic regression
var hReg = require('users/drewhart/seasonality/:fft/harmonic_regression.js');

// map stuff
var cbar = require('users/drewhart/seasonality/:viz/make_colorbar.js');
var palettes = require('users/gena/packages:palettes'); 

// QA/data-filtering assets


//--------
// IMPORTS
//--------

// read in the dataset needed
if (datasetName === 'NIRv'){
  // NIRv, its dependent variable band name, proj, scale
  var dataset = gDat.getNIRvData(minWaterOccur, nYearsNIRv);
  var dependent = 'NIRv';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
    
  // or SIF, its dependent variable band name, proj, scale
} else if (datasetName === 'SIF'){
  var dataset = gDat.getSIFData(minWaterOccur);
  var dependent = 'b1';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
  
// or TerraClimate, if calculating climatic seasonality
// (use var band name indicated by climateVar above,
// and get proj, scale, and nNeigh)
} else if (datasetName === 'TerraClimate'){
  var dataset = gDat.getTerraClimateData(minWaterOccur);
  var dependent = climateVar;
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
  
// or MODIS avg cloud proportion (for now, just averaged for all days with imagery from both Aqua and Terra)
} else if (datasetName === 'MODISCloud'){
  var dataset = gDat.getMODISCloudData(minWaterOccur);
  var dependent = 'cloud_prop';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
}



//===========
// MASKING/QA
//===========

// NOTE: only for phenology datasets
if (datasetName in {'SIF': null, 'NIRv': null}){
 
  // mask for even time series, stable LC, data missingness,
  // permutation test significance, and LC mode
    
    
  //------------------
  // filter land cover
  //------------------
  // NOTE: this will also update the mask use the preLCMask loaded above
  // NOTE: this will also reduce the temporal extent of the data
  //       if the data's original temporal extent falls outside that of the LC data!
  //       (so that we're not calling updateMask multiple separate times, since profiling
  //        showed that it chewed up a significant percent of total runtime)
  var datasetFullyMasked = maskLC.maskLandCover(maskingMode,
                                                dataset,
                                                minPctValidLC,
                                                proj,
                                                preLCMask);
                                                
  } else {
  // NOTE: skip masking, if processing climate datasets
  var datasetFullyMasked = dataset;
}


//----------------------
// CALCULATE SEASONALITY
//----------------------

// get a one-year daily time series of fitted values, derived from that regression
var maskedReg = ee.ImageCollection(hReg.calcHarmonicRegression(datasetFullyMasked, 
                                                               dependent,
                                                               harmonics, 
                                                               bandsForViz));



//----------------------------
// PRINT AND MAP INTERMEDIATES
//----------------------------

// print, if needed
if (verbose){
  print('VERBOSE OUTPUT ABOUT DATASET: ');
  print(datasetName, dataset.first());
  print('projection: ', proj);
  print('scale: ', scale);
  print('DatasetFullyMasked: ', datasetFullyMasked);
  print('------------------------------------');
}

// map intermediates, if needed
if (map_intermediates){
  Map.addLayer(dataset, {}, 'unmasked');
  Map.addLayer(datasetLooseMasked, {}, 'LooseMasked');
  Map.addLayer(datasetLCMasked, {}, 'LCMasked');
  Map.addLayer(datasetLCTSMasked, {}, 'LCTSMasked');
  Map.addLayer(datasetLCTSSigMasked, {}, 'LCTSSigMasked');
}


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
  
var output = maskedReg
  .first()
  .select(['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2', 'R2']);

var scale = dataset
  .first()
  .projection()
  .nominalScale()
  .getInfo();

// calculate the needed kernelSize in order to get adequate patch-overlap margins that allow for
// calculation of asynchrony all the way out to the neighborhood-radius distance
// NOTE: multiply by 2 because overlap margin is of size kernelSize/2
// DETH: 10-09-21: JUST REALIZED! this is not providing adequate marginal overlap
//                 at higher latitudes because I'm using the nominal scale in meters
//                 **AT THE POINT OF TRUE SCALE**, which is near the equator, not toward the poles!
//                 changing this instead to hard-code for the minimum per-cell distance in our dataset
//                 (i.e. the longitudinal distance across a cell all the way up/down at 60 N/S,
//                 given that we're working in a fixed res of ~0.05 degrees for most data,
//                 and a fixed res of ~0.04166666... for the TerraClimate data),
//                 which according to https://en.wikipedia.org/wiki/Longitude
//                 and confirmed in other sources
//                 is ~55.8 km per degree * 0.05 degrees = ~2.79km per cell for most data,
//                 and ~55.8 * 0.041666666666666664 = ~2.325km per cell for TerraClimate data,
//                 then rounding down further to 2.25km, to ensure that we take more than adequate margin sizes...
var kernelSize = maxNeighDist
  //.divide(scale)
  .divide(ee.Number(2.25*1000))
  .multiply(kernelSizeBufferFactor)
  .multiply(2)
  .ceil()
  .getInfo();
// add 1 to the kernelSize if it's odd (since it needs to be evenly divisible by 2)
if (kernelSize % 2 == 1){
  var kernelSize = ee.Number(kernelSize).add(1).getInfo();
}

if (verbose){
  print('VERBOSE OUTPUT ABOUT EXPORT PARAMETERS: ');
  print('masking mode: ', maskingMode);
  print('kernelSize', kernelSize);
  print('scale', scale);
  print('maxNeighDist', maxNeighDist);
  if (fileFormat in {'TFRecord': null, 'GeoTIFF': null}){
    print('saving to BDrive folder: ', folder);
  }
  print('------------------------------------');
}

if (datasetName != 'TerraClimate'){
  var taskAndAssetName = regionAbbrev.cat('_seas_coeffs_').cat(datasetName).getInfo();      
} else {
  var taskAndAssetName = regionAbbrev.cat('_seas_coeffs_').cat(datasetName).cat('_').cat(climateVar).getInfo();
}
if (maskingMode == 'strict'){
  var taskAndAssetName = ee.String(taskAndAssetName).cat('_STRICTMASK').getInfo();
} else if (maskingMode == 'default'){
  var taskAndAssetName = taskAndAssetName;
}
 

// Export the image, specifying scale and region.
if (fileFormat == 'TFRecord'){
  Export.image.toDrive({
    image: output,
    description: taskAndAssetName,
    folder: folder,
    scale: dataset.first().projection().nominalScale().getInfo(),
    region: roi,
    // use the file format provided at the top of the script
    fileFormat: fileFormat,
    maxPixels: maxPixels,
    formatOptions: {
      maxFileSize: maxFileSize,
      defaultValue: -9999,
      // NOTE: It appears that both patchDimensions AND kernelSize need to be expressed in
      // pixel count, not in meters...
      patchDimensions: [patchSize, patchSize],
      // NOTE: kernelSize should be such that kernelSize/2 > maxDist used in asynch calc;
      //       because in pixel count, not meters, if I need >300km overlap around edge
      //       (because that's my target maxDist that defines the neighborhood within which
      //       asynchrony is calculated), that's >300,000m = 1/2 kernel Size,
      //       so 600,000m = kernel size in meters,
      //       so 600,000/6000 (nominal scale in meters, rounded up from 5565.974539664152m)
      //       = 100 pixels, so that's my kernelSize that will ensure an adequate margin to
      //       calculate asynchrony from each patch without crosstalk between patches
      kernelSize: [kernelSize, kernelSize]
     }
   });
}
else if (fileFormat == 'Asset'){
    Export.image.toAsset({
    image: output,
    description: taskAndAssetName,
    assetId: taskAndAssetName,
    scale: dataset.first().projection().nominalScale().getInfo(),
    region: roi,
    });
}
else if (fileFormat == 'GeoTIFF'){
  Export.image.toDrive({
    image: output,
    description: taskAndAssetName,
    folder: folder,
    scale: dataset.first().projection().nominalScale().getInfo(),
    region: roi,
    // use the file format provided at the top of the script
    fileFormat: fileFormat,
    maxPixels: maxPixels,
    formatOptions: {
      defaultValue: -9999,
      // NOTE: It appears that both patchDimensions AND kernelSize need to be expressed in
      // pixel count, not in meters...
      patchDimensions: [patchSize, patchSize],
      // NOTE: kernelSize should be such that kernelSize/2 > maxDist used in asynch calc;
      //       because in pixel count, not meters, if I need >300km overlap around edge
      //       (because that's my target maxDist that defines the neighborhood within which
      //       asynchrony is calculated), that's >300,000m = 1/2 kernel Size,
      //       so 600,000m = kernel size in meters,
      //       so 600,000/6000 (nominal scale in meters, rounded up from 5565.974539664152m)
      //       = 100 pixels, so that's my kernelSize that will ensure an adequate margin to
      //       calculate asynchrony from each patch without crosstalk between patches
      kernelSize: [kernelSize, kernelSize]
    }
  });
}

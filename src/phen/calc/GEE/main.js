// This is the main script used to preprocess data for our study of
// global land surface phenology and phenological asynchrony patterns
// (Terasaki Hart et al. 2024)


//------------
// LOAD PARAMS
//------------
var params = require('users/drewhart/seasonality/:params.js');


//-----------------------------------
// LOAD REQUIRED FUNCTIONS AND ASSETS
//-----------------------------------
// data-loading
var gDat = require('users/drewhart/seasonality/:io/get_data.js');

// harmonic regression
var hReg = require('users/drewhart/seasonality/:fft/run_harmonic_regression.js');


//--------
// IMPORTS
//--------
// read in the dataset needed
if (params.datasetName === 'NIRv'){
  // NIRv, its dependent variable band name, proj, scale
  var dataset = gDat.getNIRvData(params.maskLowBRDFAlbQual,
                                 params.maxBRDFAlbQualVal,
                                 params.maskWaterSeparately,
                                 params.NIRvNYears,
                                 params.NIRvEndYear,
                                 params.NIRvDayStep,
                                 params.NIRvClampToMinPos,
                                 params.NIRvUnmaskToMinPos,
                                 params.NIRvMaskLteZero
                                );
  var dependent = 'NIRv';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
  // or SIF, its dependent variable band name, proj, scale
} else if (params.datasetName === 'SIF'){
  var dataset = gDat.getSIFData(params.maskWaterSeparately,
  // NOTE: match min-positive clamping and unmasking approach used for NIRv data
                                params.NIRvClampToMinPos,
                                params.NIRvUnmaskToMinPos
                               );
  var dependent = 'b1';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
  
// or TerraClimate, if calculating climatic seasonality
// (use var band name indicated by climateVar above,
// and get proj, scale, and nNeigh)
} else if (params.datasetName === 'TerraClimate'){
  var dataset = gDat.getTerraClimateData(params.maskWaterSeparately);
  var dependent = params.climateVar;
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
  
// or MODIS avg cloud proportion (for now, just averaged for all days with imagery from both Aqua and Terra)
} else if (params.datasetName === 'MODISCloud'){
  var dataset = gDat.getMODISCloudData();
  var dependent = 'cloud_prop';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
}


//---------------
// GET REGRESSION
//---------------
if (params.datasetName != 'NIRv'){
  // get a one-year daily time series of fitted values, derived from that regression
  var reg = ee.ImageCollection(hReg.calcHarmonicRegression(dataset, 
                                                           dependent,
                                                           params.harmonics, 
                                                           params.bandsForViz));
  // grab the regression output's first Image, which is the result (has coefficients and R2)
  var result = reg
    .first()
    .select(['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2', 'R2']);
} else {
  var result = ee.Image('users/drewhart/NIRv_harmonic_regression')
    .select(['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2', 'R2']);
}


//----------------
// MASK THE RESULT
//----------------
// NOTE: we only mask the phenology datasets
if (params.datasetName in {'SIF': null, 'NIRv': null}){
  // load the mask
  if (params.maskingMode == 'default'){
    var mask = ee.Image('users/drewhart/LSP_mask_DEFAULT');
  } else if (params.maskingMode == 'strict'){
    var mask = ee.Image('users/drewhart/LSP_mask_STRICT');
  }
  // apply it
  var result = result.updateMask(mask);
}




//----------------------------
// PRINT AND MAP INTERMEDIATES
//----------------------------
// print, if needed
if (params.verbose){
  print('VERBOSE OUTPUT ABOUT DATASET: ');
  print('size of dataset: ', dataset.size());
  print(params.datasetName, dataset.first());
  print('result', result);
  print('------------------------------------');
}
if (params.map_intermediates){
  Map.setOptions('SATELLITE');
  Map.addLayer(dataset, {}, 'dataset');
  if (params.datasetName in {'SIF': null, 'NIRv': null}){
    Map.addLayer(mask, {min:0, max:1, palette:['red', 'white']}, 'mask');
  }
  if (params.datasetName != 'NIRv'){
    Map.addLayer(ee.ImageCollection(reg.toList(reg.size().getInfo()-1, 1)).select('fitted'), {}, 'reg');
  }
  Map.addLayer(result.select('R2'), {}, 'R2');
}




//----------
// EXPORT IT
//----------
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
//                 (i.e. the longitudinal distance across a cell all the way up/down,
//                 from 60N to 60 S (NOTE: we actually go to 75N but data drops out just before 60N
//                 because of the avg monthly proportion and monthly evenness masks,
//                 given that we're working in a fixed res of ~0.05 degrees for most data,
//                 and a fixed res of ~0.04166666... for the TerraClimate data),
//                 which according to https://en.wikipedia.org/wiki/Longitude
//                 and confirmed in other sources
//                 is ~55.8 km per degree * 0.05 degrees = ~2.79km per cell for most data,
//                 and ~55.8 * 0.041666666666666664 = ~2.325km per cell for TerraClimate data,
//                 then rounding down further to 2.25km, to ensure that we take more than adequate margin sizes...
var kernelSize = params.maxNeighDist
  //.divide(scale)
  .divide(ee.Number(2.25*1000))
  .multiply(params.kernelSizeBufferFactor)
  .multiply(2)
  .ceil()
  .getInfo();
// add 1 to the kernelSize if it's odd (since it needs to be evenly divisible by 2)
if (kernelSize % 2 == 1){
  var kernelSize = ee.Number(kernelSize).add(1).getInfo();
}

if (params.verbose){
  print('VERBOSE OUTPUT ABOUT EXPORT PARAMETERS: ');
  print('masking mode: ', params.maskingMode);
  print('patchDimensions: ', [params.patchSize, params.patchSize]);
  print('kernelSize', kernelSize);
  print('projection: ', proj);
  print('scale', scale);
  print('maxNeighDist', params.maxNeighDist);
  if (params.fileFormat in {'TFRecord': null, 'GeoTIFF': null}){
    print('saving to BDrive folder: ', params.mainFolder);
  }
  print('------------------------------------');
}

if (params.datasetName != 'TerraClimate'){
  var taskAndAssetName = params.regionAbbrev.cat('_coeffs_').cat(params.datasetName).getInfo();      
} else {
  var taskAndAssetName = params.regionAbbrev.cat('_coeffs_').cat(params.datasetName).cat('_').cat(params.climateVar).getInfo();
}
if (params.maskingMode == 'strict'){
  var taskAndAssetName = ee.String(taskAndAssetName).cat('_STRICTMASK').getInfo();
} else if (params.maskingMode == 'default'){
  var taskAndAssetName = taskAndAssetName;
}
 

// Export the image, specifying scale and region.
if (params.fileFormat == 'TFRecord'){
  Export.image.toDrive({
    image: result,
    description: taskAndAssetName,
    folder: params.mainFolder,
    scale: dataset.first().projection().nominalScale().getInfo(),
    region: params.roi,
    fileFormat: params.fileFormat,
    maxPixels: params.maxPixels,
    formatOptions: {
      maxFileSize: params.maxFileSize,
      defaultValue: params.defaultValue,
      // NOTE: both patchDimensions and kernelSize are expressed in pixel count, not in meters
      patchDimensions: [params.patchSize, params.patchSize],
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
else if (params.fileFormat == 'Asset'){
    Export.image.toAsset({
    image: result,
    description: taskAndAssetName,
    assetId: taskAndAssetName,
    scale: dataset.first().projection().nominalScale().getInfo(),
    region: params.roi,
    });
}
else if (params.fileFormat == 'GeoTIFF'){
  Export.image.toDrive({
    image: result,
    description: taskAndAssetName,
    folder: params.mainFolder,
    scale: dataset.first().projection().nominalScale().getInfo(),
    region: params.roi,
    fileFormat: params.fileFormat,
    maxPixels: params.maxPixels,
  });
}


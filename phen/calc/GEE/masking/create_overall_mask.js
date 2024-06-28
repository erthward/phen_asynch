// Make mask map for both LSP datasets and for both for default (all vals <2) and strict (all vals <1) masking modes,
// then export assets for both of those modes

//=======
// PARAMS
//=======

var params = require('users/drewhart/seasonality/:params.js');


//==========================
// LOAD FUNCTIONS AND ASSETS
//==========================

var lcMask = ee.Image('users/drewhart/LSP_lcMask');
var tsPctDataAvailability = ee.Image('users/drewhart/NIRv_ts_pct_data_availability');
var tsEvenness = ee.Image('users/drewhart/NIRv_monthly_Pielou_evenness')
  // NOTE: some pixels are getting values erroneously >1, but I spent a long time inspecting and it appears
  //       only to happen offshore and I manually checked that the Pielou evenness code is working for
  //       real data at some pixels, so I feel confident in the result and will just clamp values to 1
  .clamp({low: 0, high: 1});
var signifMask = ee.Image('users/drewhart/NIRv_permutation_significance_mask');
var gDat = require('users/drewhart/seasonality/:io/get_data.js');




//============
// DATA IMPORT
//============

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
  var dataset = gDat.getSIFData(params.maskWaterSeparately);
  var dependent = 'b1';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
}




//===========
// MASKING/QA
//===========

//--------------------------
// impose evenness threshold
//--------------------------
// NOTE: unmask to explicitly include zeros in the overallMask in places where the
// min avg monthly proportional data availability mask actually masked pixels
// out of the evenness mask
var evennessMask = tsEvenness.unmask().gte(params.minTSEvenness);


//--------------------
// impose LC threshold
//--------------------
if (params.maskingMode == 'default'){
  // get all pixels that are valid (i.e., not barren, ice, water, urban)
  // and either always 'natural' or always ag across our time series...
  var lcMaskThresholded = lcMask.gt(0);
} else {
  // ... or get only valid pixels that are always 'natural' across our time series
  var lcMaskThresholded = lcMask.gt(1);
}
  
  
//-----------------------------
// impose missingness threshold
//-----------------------------
var shortTSMask = tsPctDataAvailability.gte(params.minPctTSDataAvailability);


//---------------------
// compose overall mask
//---------------------
var overallMask = evennessMask.and(lcMaskThresholded).and(shortTSMask).and(signifMask);
if (params.map_intermediates){
  Map.addLayer(overallMask, {min:0, max:1, palette:['red', 'white']}, 'overallMask');
}




//====================
// EXPORT OVERALL MASK
//====================
var scale = dataset
  .first()
  .projection()
  .nominalScale()
  .getInfo();

if (params.maskingMode == 'strict'){
  var taskAndAssetName = 'LSP_mask_STRICT';
} else if (params.maskingMode == 'default'){
  var taskAndAssetName = 'LSP_mask_DEFAULT';
}
  
// Export the image, specifying scale and region.

Export.image.toAsset({
image: overallMask,
description: taskAndAssetName,
assetId: taskAndAssetName,
scale: scale,
region: params.roi,
});




//=============================================
// PLOT AND EXPORT INDIVIDUAL MASKS AS GEOTIFFS
//=============================================

// NOTE: month_props_min_mask separately exported by calc_MODIS_month_evenness.js

if (params.plotAndExportMasks && params.datasetName == 'NIRv'){
  // grab the water mask,
  // export that to asset,
  // and also mask other masks using the water mask and then export them
  var waterMask = ee.Image('users/drewhart/LSP_waterMask');
  Map.addLayer(evennessMask.unmask().and(waterMask).selfMask(), {palette:['blue']}, 'evennessMask', false);
  Map.addLayer(lcMask.unmask().and(waterMask).selfMask(), {palette:['black', 'yellow', 'red']}, 'lcMask (default and strict)', false);
  Map.addLayer(shortTSMask.unmask().and(waterMask).selfMask(), {palette:['blue']}, 'shortTSMask', false);
  Map.addLayer(signifMask.unmask().and(waterMask).selfMask(), {palette:['blue']}, 'signifMask', false);
  Export.image.toDrive({image: waterMask,
                        description: 'waterMask',
                        folder: 'LSP_mask_outputs_from_GEE',
                        region: params.roi,
                        scale: scale,
                        maxPixels: params.maxPixels,
                        shardSize: params.shardSize,
                        fileDimensions: params.fileDimensions,
                        fileFormat: 'GeoTIFF'
                        }
                       );
  Export.image.toDrive({image: lcMask.unmask().updateMask(waterMask),
                       description: 'lcMask',
                       folder: 'LSP_mask_outputs_from_GEE',
                       region: params.roi,
                       scale: scale,
                       maxPixels: params.maxPixels,
                       shardSize: params.shardSize,
                       fileDimensions: params.fileDimensions,
                       fileFormat: 'GeoTIFF'
                       }
                      );
  Export.image.toDrive({image: shortTSMask.unmask().updateMask(waterMask),
                       description: ee.String('shortTSMask_')
                                      .cat(ee.String(params.datasetName)).getInfo(),
                       folder: 'LSP_mask_outputs_from_GEE',
                       region: params.roi,
                       scale: scale,
                       maxPixels: params.maxPixels,
                       shardSize: params.shardSize,
                       fileDimensions: params.fileDimensions,
                       fileFormat: 'GeoTIFF'
                       }
                      );
  Export.image.toDrive({image: evennessMask.unmask().updateMask(waterMask),
                       description: ee.String('evennessMask_')
                                      .cat(ee.String(params.datasetName)).getInfo(),
                       folder: 'LSP_mask_outputs_from_GEE',
                       region: params.roi,
                       scale: scale,
                       maxPixels: params.maxPixels,
                       shardSize: params.shardSize,
                       fileDimensions: params.fileDimensions,
                       fileFormat: 'GeoTIFF'
                       }
                      );
  Export.image.toDrive({image: signifMask.unmask().updateMask(waterMask),
                       description: ee.String('signifMask_')
                                      .cat(ee.String(params.datasetName)).getInfo(),
                       folder: 'LSP_mask_outputs_from_GEE',
                       region: params.roi,
                       scale: scale,
                       maxPixels: params.maxPixels,
                       shardSize: params.shardSize,
                       fileDimensions: params.fileDimensions,
                       fileFormat: 'GeoTIFF'
                       }
                      );                     
}
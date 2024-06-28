// params and imports
var params = require('users/drewhart/seasonality/:params.js');
var gDat = require('users/drewhart/seasonality/:io/get_data.js');

// create an ImageCollection from all 12 monthly maps of average proportional data-availability
var month_props = ee.ImageCollection(['users/drewhart/NIRv_data_avail_prop_month01',
                                      'users/drewhart/NIRv_data_avail_prop_month02',
                                      'users/drewhart/NIRv_data_avail_prop_month03',
                                      'users/drewhart/NIRv_data_avail_prop_month04',
                                      'users/drewhart/NIRv_data_avail_prop_month05',
                                      'users/drewhart/NIRv_data_avail_prop_month06',
                                      'users/drewhart/NIRv_data_avail_prop_month07',
                                      'users/drewhart/NIRv_data_avail_prop_month08',
                                      'users/drewhart/NIRv_data_avail_prop_month09',
                                      'users/drewhart/NIRv_data_avail_prop_month10',
                                      'users/drewhart/NIRv_data_avail_prop_month11',
                                      'users/drewhart/NIRv_data_avail_prop_month12'
                                     ]);
                   
// mask out pixels for which any month's mean prop of valid days is below the min
var month_props_min_mask = month_props
// NOTE: first unmask all Images, to avoid the banding that otherwise occurs across
// northern latitudes because months with average proportional data availability equal to 0.0
// get masked out and fall out of minimum threshold assessment across all months
  .map(function(img){return img.unmask()})
  .reduce(ee.Reducer.min())
  .gte(params.minMonthlyPropDataAvail);
var month_props_min_masked = month_props.map(function(img){
  return img.updateMask(month_props_min_mask)});
  
//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
// export the minimum average monthly proportion data-availability mask Image
var waterMask = gDat.maskWaterMODIS(ee.ImageCollection(ee.Image.constant(1)), 100).first();
var scale = gDat.getSIFData(0.5).first().projection().nominalScale().getInfo();
if (params.map_intermediates){
  Map.addLayer(month_props, {min: 0, max:1 }, 'month_props', false);
  Map.addLayer(month_props_min_mask, {min: 0, max:1 }, 'month_props_min_mask');
}
Export.image.toDrive({image: month_props_min_mask.and(waterMask),
                      description: 'monthPropsMinMask_NIRv',
                      folder: 'LSP_mask_outputs_from_GEE',
                      region: params.roi,
                      scale: scale,
                      maxPixels: params.maxPixels,
                      shardSize: params.shardSize,
                      fileDimensions: params.fileDimensions,
                      fileFormat: 'GeoTIFF'
                      }
                      );
//. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

// Image: sum of all months' average proportions
var month_props_sum = month_props_min_masked.reduce(ee.Reducer.sum());
if (params.map_intermediates){
  Map.addLayer(month_props_sum, {}, 'month_props_sum');
}

// ImageCollection: proportion of each pixel's annual time series that is represented by each month
var month_ann_props = month_props_min_masked
  .map(function(img){return img.divide(month_props_sum).rename('prop')});

// only keep pixels that are unmasked (i.e., have at least some data) for all 12 months
var month_ann_props_valid = month_ann_props.reduce(ee.Reducer.count()).eq(12);
var month_ann_props = month_ann_props.map(function(img){return img.updateMask(month_ann_props_valid)});
if (params.map_intermediates){
  Map.addLayer(month_ann_props, {}, 'month_ann_props');
}

// Image: Shannon diversity of each pixel's time series
var shannon_div = month_ann_props
  .map(function(img){return img.multiply(img.log())})
  .reduce(ee.Reducer.sum())
  .multiply(ee.Image.constant(-1));
// Image: max possible Shannon diversity (constant: ln(12))
var max_shannon_div = ee.Image.constant(ee.Number(12).log());
// Image: Pielou's evenness of each pixel's time series
// (NOTE: calculated by months, so not using perfectly even divisions of the calendar year, but close enough!)
var month_evenness = shannon_div.divide(max_shannon_div).toFloat();
if (params.map_intermediates){
  Map.addLayer(month_evenness, {min:0, max:1}, "Pielou's evenness of NIRv time series");
}

// export as a masking layer
Export.image.toAsset({image: month_evenness,
                      description: 'NIRv_monthly_Pielou_evenness',
                      assetId: 'users/drewhart/NIRv_monthly_Pielou_evenness',
                      region: params.roi,
                      pyramidingPolicy: 'mean',
                      maxPixels: params.maxPixels,
                      scale: month_props.first().projection().nominalScale().getInfo()});

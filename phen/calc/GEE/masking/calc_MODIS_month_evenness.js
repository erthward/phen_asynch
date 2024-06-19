// minimum mean monthly proportion of valid days
// (below which to mask out a pixel for not having adequate monthly coverage)
var minMonthlyPropVal = 0.1;

var gDat = require('users/drewhart/seasonality/:io/get_data.js');
// NOTE: 8 YEARS IS THE MAX LENGTH GEE WOULD CALCULATE, BUT THAT SHOULD BE FINE
//       BECAUSE IT MOST LIKELY ONLY MAKES THE END RESULT
//       MORE CONSERVATIVE THAN IF WE HAD 10 YEARS
var NIRv = gDat.getNIRvData(50, 10)
    .filterDate('2010-01-01T00:00', '2019-12-31T00:00');
print('NIRv time series', NIRv);
// var num_im_to_map = 1;
// var trueColor = ee.Image(ee.List(MODIS.toList(1, num_im_to_map)).get(0)).select([
//   'Nadir_Reflectance_Band1', 'Nadir_Reflectance_Band4',
//   'Nadir_Reflectance_Band3'
// ]);
// var trueColorVis = {
//   min: 0.0,
//   max: 4000.0,
//   gamma: 1.4,
// };
// Map.addLayer(trueColor, trueColorVis, 'first MODIS image (true color)');

// calc average number of values in each month
var months = ee.List(['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']);
var month_lengths = ee.Dictionary({'01': 31,
                 '02': 28,
                 '03': 31,
                 '04': 30,
                 '05': 31,
                 '06': 30,
                 '07': 31,
                 '08': 31,
                 '09': 30,
                 '10': 31,
                 '11': 30,
                 '12': 31});
// ImageCollection: average proportion of each month's days having valid data in a pixel's time series                 
var month_props = ee.ImageCollection(months
  .map(function(month){return NIRv
    .map(function(img){return img
      .select('NIRv')
      .set({'month': ee.Number.parse(ee.String(img.get('system:index')).slice(5,7))})})
    .filter(ee.Filter.expression(ee.String('month==').cat(month)))
    .reduce(ee.Reducer.count())
    .divide(ee.Image.constant(10))
    .divide(ee.Image.constant(ee.Dictionary(month_lengths).get(month)))
    // NOTE: clamp to a minimum of 0.0001 (<<1/31), to avoid division-by-zero issues downstream
    .clamp(0.0001, 1);
    }));
Map.addLayer(month_props, {}, 'month_props');

// mask out pixels for which any month's mean prop of valid days is below the min
var month_props_min_mask = month_props
  .reduce(ee.Reducer.min())
  .gte(minMonthlyPropVal);
var month_props_min_masked = month_props.map(function(img){
  return img.updateMask(month_props_min_mask)});

// Image: sum of all months' average proportions
var month_props_sum = month_props_min_masked.reduce(ee.Reducer.sum());
Map.addLayer(month_props_sum, {}, 'month_props_sum');

// ImageCollection: proportion of each pixel's annual time series that is represented by each month
var month_ann_props = month_props_min_masked
  .map(function(img){return img.divide(month_props_sum)});

// only keep pixels that are unmasked (i.e., have at least some data) for all 12 months
var month_ann_props_valid = month_ann_props.reduce(ee.Reducer.count()).eq(12);
var month_ann_props = month_ann_props.map(function(img){return img.updateMask(month_ann_props_valid)});
Map.addLayer(month_ann_props, {}, 'month_ann_props');

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
Map.addLayer(month_evenness, {min:0, max:1}, "Pielou's evenness of NIRv time series");

// export as a masking layer
Map.setZoom(0); // export globally
//var roi = ee.Geometry.Polygon(
//        [[[-76.93671917189324, 6.885058085787389],
//          [-76.93671917189324, 2.8136426929061034],
//          [-71.24580120314324, 2.8136426929061034],
//          [-71.24580120314324, 6.885058085787389]]], null, false);
Export.image.toAsset({image: month_evenness,
                      description: 'NIRv_monthly_Pielou_evenness_10YR',
                      assetId: 'users/drewhart/NIRv_monthly_Pielou_evenness_10YR',
                      //region: roi,
                      pyramidingPolicy: 'mean',
                      maxPixels: 4000000000,
                      scale: NIRv   .first().projection().nominalScale().getInfo()});

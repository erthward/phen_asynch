// set NIRv-calculating params and read NIRv dataset in
var params = require('users/drewhart/seasonality/:params.js');
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
var NIRv = gDat.getNIRvData(params.maskLowBRDFAlbQual,
                            params.maxBRDFAlbQualVal,
                            params.maskWaterSeparately,
                            params.NIRvNYears,
                            params.NIRvEndYear,
                            params.NIRvDayStep,
                            false,  // do not need to bother clamping to the minimum positive value observed
                            params.NIRvUnmaskToMinPos,
                            false   // do not want to mask values <= 0, since we're not clamping to min pos first
                           );
print('NIRv time series', NIRv);

// add a month property to all Images in the NIRv ImageCollection
var NIRvWithMonth = NIRv
  .map(function(img){return img
    .set({'month': ee.Number.parse(ee.String(img.get('system:index')).slice(5,7))});
  });
if (false){
//if (params.map_intermediates){
  Map.addLayer(NIRvWithMonth, {}, 'NIRvWithMonth', false);
}

// calc average number of values in each month
var months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'];
var month_lengths = ee.Dictionary({'01': 31/params.NIRvDayStep,
                 '02': 28/params.NIRvDayStep,
                 '03': 31/params.NIRvDayStep,
                 '04': 30/params.NIRvDayStep,
                 '05': 31/params.NIRvDayStep,
                 '06': 30/params.NIRvDayStep,
                 '07': 31/params.NIRvDayStep,
                 '08': 31/params.NIRvDayStep,
                 '09': 30/params.NIRvDayStep,
                 '10': 31/params.NIRvDayStep,
                 '11': 30/params.NIRvDayStep,
                 '12': 31/params.NIRvDayStep});

// for each month, create an asset reporting the average proportion
// of each month's days in a pixel's time series for which data is available
for (var m = 0; m < months.length; m++){
  var month = months[m];
  var month_prop = NIRvWithMonth
      .filter(ee.Filter.expression(ee.String('month==').cat(ee.String(month))))
      // NOTE: updateMask(1) will change all partially transparent pixels
      //       (i.e., pixels with <100% available native-res data within them) to opaque(
      //       (by setting all non-zero masks values to 1)
      .map(function(img){return img.multiply(0).add(1).updateMask(1)})
      .sum()
      .divide(ee.Image.constant(params.NIRvNYears).toFloat())
      .divide(ee.Image.constant(ee.Dictionary(month_lengths).get(ee.String(month))).toFloat())
      // NOTE: clamp to a minimum of 0.0001 (<<1/31), to avoid division-by-zero issues downstream
      .clamp(0.0001, 1);
  if (params.map_intermediates && m===0){
    Map.addLayer(NIRvWithMonth.filter(ee.Filter.expression(ee.String('month==').cat(ee.String(month)))), {}, 'month-filtered NIRv', false);
    Map.addLayer(NIRvWithMonth.filter(ee.Filter.expression(ee.String('month==').cat(ee.String(month)))).map(function(img){return img.multiply(0).add(1)}), {}, 'month-filtered NIRv masks', false);
    Map.addLayer(month_prop, {}, ee.String('month_prop: ').cat(ee.String(month)).getInfo(), false);
  }
  // export to an asset
  // export as a masking layer
  var assetId = ee.String('NIRv_data_avail_prop_month').cat(ee.String(month));
  Export.image.toAsset({image: month_prop,
                        description: assetId.getInfo(),
                        assetId: ee.String('users/drewhart/').cat(assetId).getInfo(),
                        region: params.roi,
                        pyramidingPolicy: 'mean',
                        maxPixels: params.maxPixels,
                        scale: NIRv.first().projection().nominalScale().getInfo()
                      });
}
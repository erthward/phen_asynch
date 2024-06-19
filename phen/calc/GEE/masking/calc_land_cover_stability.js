var lc = ee.ImageCollection("MODIS/006/MCD12Q1")
  .select('LC_Type1')
  .filterDate('2009-01-01T00:00', '2019-12-31T00:00');
print(lc);

// determine mode of each pixel's LC type across the 10-year time series
// NOTE: using mode to avoid >1 mode
// NOTE: `maxRaw: 20` ensures no histogram binning used to calculate mode, because 20 > 10 years and > 17 LC types
var majority = lc.reduce(ee.Reducer.mode({maxRaw: 20}));

// determine percent of time series each pixel spends in majority LC
var pct_time_lc_mode = lc
  .map(function(img){return img.eq(majority)})
  .reduce(ee.Reducer.sum())
  .divide(lc.reduce(ee.Reducer.count()));
Map.addLayer(pct_time_lc_mode, {min: 0, max: 1}, 'percent time in LC mode');

// export as a masking layer
Map.setZoom(0); // export globally
Export.image.toAsset({image: pct_time_lc_mode,
                      description: 'pct_time_in_MODIS_lc_mode',
                      assetId: 'users/drewhart/pct_time_in_MODIS_lc_mode',
                      pyramidingPolicy: 'mean',
                      maxPixels: 4000000000,
                      scale: lc.first().projection().nominalScale().getInfo()});

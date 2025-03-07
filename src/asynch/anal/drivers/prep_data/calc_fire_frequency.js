// get params
var params = require('users/drewhart/seasonality/:params.js');

// get SIF dataset and projection
var SIF = ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN');
var proj = SIF.first().projection();

// load MODIS monthly burned area dataset
var burnedArea = ee.ImageCollection('MODIS/061/MCD64A1');
print(burnedArea);

// use burndate to count unique burn events
var origProj = burnedArea.first().projection();
var burnCount = burnedArea
   .select('BurnDate')
   .reduce(ee.Reducer.count())
   // NOTE: projection gets dropped by reduce, so needs to be reassigned
   .setDefaultProjection(origProj);
   
// unmask (retaining zeros where burns are unobserved),
// then divide by total length (max potential number of
// burn events identified) to express as burn frequency
var burnFreq = burnCount.unmask().divide(burnedArea.size());
if (params.map_intermediates){
  Map.addLayer(burnFreq, {palette: ['black', 'red', 'orange', 'yellow', 'white'], min:0, max:1/25, opacity:0.8}, 'burnFreq');
}

// reduce to match SIF res
var burnFreq_agg = burnFreq
  .reduceResolution(ee.Reducer.mean(), false, 250)
  .reproject(proj);
if (params.map_intermediates){
  Map.addLayer(burnFreq_agg, {palette: ['black', 'red', 'orange', 'yellow', 'white'], min:0, max:1/25, opacity:0.8}, 'burnFreq_agg'); 
}

// mask to only pixels whose LSP signal is being factored into the LSP asynchrony map
var strict_mask = ee.Image('users/drewhart/LSP_mask_STRICT');
var burnFreq_agg_mask = burnFreq_agg.mask(strict_mask);
if (params.map_intermediates){
  Map.addLayer(burnFreq_agg_mask, {palette: ['black', 'red', 'orange', 'yellow', 'white'], min:0, max:1/25, opacity:0.8}, 'burnFreq_agg_mask'); 
}

// calc neighborhood mean within 100km
var burnFreq_mean = burnFreq_agg_mask
  .reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.circle(100000, 'meters'));
if (params.map_intermediates){
  Map.addLayer(burnFreq_mean, {palette: ['black', 'red', 'orange', 'yellow', 'white'], min:0, max:1/25, opacity:0.8}, 'burnFreq_mean');
}

// export
var scale = proj
  .nominalScale()
  .getInfo();
print('scale', scale);

Export.image.toDrive({
  image: burnFreq_mean,
  description: 'MODIS_fire_freq_mean',
  folder: 'LSP_asynch_drivers_analysis_covars',
  scale: scale,
  region: params.roi,
  fileFormat: 'GeoTIFF',
  maxPixels: params.maxPixels,
  shardSize: params.shardSize,
  fileDimensions: params.fileDimensions,
});

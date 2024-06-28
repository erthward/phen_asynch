// get params
var params = require('users/drewhart/seasonality/:params.js');

// get SIF dataset and projection
var SIF = ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN');
var proj = SIF.first().projection();

// load MODIS yearly LC data, and
// filter to same years as NIRv dataset
// NOTE: switched from MCD12Q1 to MCD12C1, per Reviewer 1's recommendation
var lc = ee.ImageCollection("MODIS/061/MCD12C1")
  .select('Majority_Land_Cover_Type_1')
  .filterDate(ee.String(ee.Number(params.NIRvEndYear - params.NIRvNYears - 1)).cat('-12-31T00:00'),
              ee.String(ee.Number(params.NIRvEndYear -1)).cat('-12-31T00:00'));
if (params.verbose){
  print('lc', lc);
}

// reclass to reflect vegetation structure
var lc_reclass = lc.map(
  function(img){return img.multiply(0)
    .where(img.gt(0).and(img.lt(6)), 1)     // forest = 1
    .where(img.eq(6).or(img.eq(7)), 2)      // shrubland = 2
    .where(img.eq(8).or(img.eq(9)), 3)      // savanna/woodland = 3
    .where(img.eq(10), 4)                   // grassland = 4
    .where(img.eq(11), 5)                   // perm. wetland = 5
    // NOTE: same as in strict LC masking, mask out ag, urban, barren, ice
    // NOTE: also permanent water, which is already zero
    .where(img.gte(12), 0)                  
    .selfMask();
  });
if (params.map_intermediates){
  var lc_viz = {min:1, max:5, palette: ['darkgreen', 'brown', 'yellow', 'green', 'blue']};
  Map.addLayer(lc_reclass, lc_viz, 'lc_reclass', false);
}

// reduce to a single image using mode,
// then take just the IGBP dataset
var lc_img = lc_reclass
  .reduce(ee.Reducer.mode())
  .select('Majority_Land_Cover_Type_1_mode')
  .rename('IGBP');
if (params.map_intermediates){
  Map.addLayer(lc_img, lc_viz, 'lc_img'); 
}

// reduce to SIF res
var lc_img_agg = lc_img
  .setDefaultProjection(lc.first().projection())
  .reduceResolution(ee.Reducer.mode(),false, 250)
  .reproject(proj);
if (params.map_intermediates){
  Map.addLayer(lc_img_agg, lc_viz, 'lc_img_agg'); 
}

// calc neighborhood entropy within 100km
var lc_entropy = lc_img_agg
  .entropy(ee.Kernel.circle(100000, 'meters'))
  // NOTE: clamp to [0, 1]
  .clamp({low: 0, high: 1});
if (params.map_intermediates){
  Map.addLayer(lc_entropy, {}, 'lc_entropy'); 
}


// export
var scale = proj
  .nominalScale()
  .getInfo();

print('scale', scale);

 
Export.image.toDrive({
  image: lc_entropy,
  description: 'MODIS_IGBP_veg_entropy',
  folder: 'LSP_asynch_drivers_analysis_covars',
  scale: scale,
  region: params.roi,
  fileFormat: 'GeoTIFF',
  maxPixels: params.maxPixels,
  shardSize: params.shardSize,
  fileDimensions: params.fileDimensions,
});


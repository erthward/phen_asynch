// get SIF dataset and projection
var SIF = ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN');
var proj = SIF.first().projection();

// load MODIS yearly LC data, and
// filter to same years as NIRv dataset
var lc = ee.ImageCollection("MODIS/006/MCD12Q1")
  .filterDate('2009-12-31T00:00', '2019-12-31T00:00');

// reclass to reflect vegetation structure
var lc_reclass = lc.map(
  function(img){return img.multiply(0)
    .where(img.lt(6), 1)                    // forest = 1
    .where(img.eq(6).or(img.eq(7)), 2)      // shrubland = 2
    .where(img.eq(8).or(img.eq(9)), 3)      // savanna = 3
    .where(img.eq(10), 4)                   // grassland = 4
    .where(img.eq(11), 5)                   // perm. wetland = 5
    // NOTE: same as in strict LC masking, mask out ag, urban, barren, ice, water
    .where(img.gte(12), 0)                  
    .selfMask();
  });

// reduce to a single image using mode,
// then take just the IGBP dataset
var lc_img = lc_reclass
  .reduce(ee.Reducer.mode())
  .select('LC_Type1_mode')
  .rename('IGBP');

Map.addLayer(lc_img, {}, 'lc_img'); 
  
// reduce to SIF res
var lc_img_agg = lc_img
  // I manually checked that all the original IC's Images are nearly equal (to within 4 10_000ths of a deg),
  // so just manually reset default proj to the IC's first image's default proj
  .setDefaultProjection(lc.first().projection())
  .reduceResolution(ee.Reducer.mode(),false, 250)
  .reproject(proj);

// calc neighborhood entropy within 100km
var lc_entropy = lc_img_agg
  .entropy(ee.Kernel.circle(100000, 'meters'))
  // NOTE: clamp to [0, 1]
  .clamp({low: 0, high: 1});
Map.addLayer(lc_entropy, {}, 'LC entropy');


// export
var roi =
  ee.Geometry.Polygon(
      [[[-165, 80],
        [-165, -60],
        [180, -60],
        [180, 80],
        [-165,80]]], null, false);

var scale = lc_entropy
  .projection()
  .nominalScale()
  .getInfo();

print('scale', scale);

 
Export.image.toDrive({
  image: lc_entropy,
  description: 'MODIS_IGBP_veg_entropy',
  folder: 'ancillary_outputs_from_GEE',
  scale: scale,
  region: roi,
  // use the file format provided at the top of the script
  fileFormat: 'GeoTIFF',
  maxPixels: 100000000,
  //dimensions: 7000, // max dimension is width of 7000
  //fileDimensions: [6912, 2560] // overshoots the expected 6900x2400 dims in order 
                               //to be divisble by the 256-default shard size
});

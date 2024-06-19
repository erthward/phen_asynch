// export the combined results?
var export_combo = false;

// read in the separate results
var alreadyCombined = ee.Image('users/drewhart/testRes_0-174_combo');
var res1 = ee.Image('users/drewhart/testRes_SIF_seeds_175_224');
var res2 = ee.Image('users/drewhart/testRes_SIF_seeds_225_299');

// total number of tests run thus far
var nTests = ee.Number(300);

// alpha, to determine significance (for mapping only)
var alpha = ee.Number(0.01);
var threshold = alpha.multiply(nTests).ceil().getInfo();

// add the separate results together
var combo = alreadyCombined.add(res1).add(res2);

// map it
var palettes = require('users/gena/packages:palettes'); 
Map.addLayer(combo, {min: 0, max: threshold,
                   palette: palettes.crameri.roma[10].reverse(), opacity: 0.8},
             'combined results');

// export, if needed          
if (export_combo){
  var scale = combo.projection().nominalScale();
  // NOTE: important to use the same ROI as used in creating the separate results
  var roi = 
      ee.Geometry.Polygon(
          [[[-170, 50],
            [-170, -50],
            [179, -50],
            [179, 50]]], null, false);
  
  // Export the image to an Earth Engine asset.
  var taskAndAssetName = ee.String('testRes_0-')
    .cat(ee.String(nTests.subtract(1)))
    .cat('_combo')
    .getInfo();
  Export.image.toAsset({
    image: combo,
    description: taskAndAssetName,
    assetId: taskAndAssetName,
    scale: scale.getInfo(),
    region: roi,
    maxPixels: 4e9,
    pyramidingPolicy: {
      '.default': 'mean',
    }
  });
}

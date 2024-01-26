// max p-value beyond which to mask out a pixel (i.e., alpha value)
var maxpVal = 0.01;

var dataset = 'NIRv';
//var dataset = 'SIF';

// combine all permutation fail-count maps
if (dataset == 'SIF'){
  var perms_0_24 = ee.Image('users/drewhart/SIF_permutations_seeds_0-24');
  var perms_25_49 = ee.Image('users/drewhart/SIF_permutations_seeds_25-49');
  var perms_50_74 = ee.Image('users/drewhart/SIF_permutations_seeds_50-74');
  var perms_75_99 = ee.Image('users/drewhart/SIF_permutations_seeds_75-99');
  var allPerms = perms_0_24.add(perms_25_49).add(perms_50_74).add(perms_75_99);
  var total_ct = 100;
  var output_scale = perms_0_24.projection().nominalScale();
} else {
  var perms_0_4 = ee.Image('users/drewhart/NIRv_permutations_seeds0-4');
  var perms_5_9 = ee.Image('users/drewhart/NIRv_permutations_seeds5-9');
  var perms_10_14 = ee.Image('users/drewhart/NIRv_permutations_seeds10-14');
  var perms_15_19 = ee.Image('users/drewhart/NIRv_permutations_seeds15-19');
  var perms_20_24 = ee.Image('users/drewhart/NIRv_permutations_seeds20-24');
  var allPerms = perms_0_4.add(perms_5_9).add(perms_10_14).add(perms_15_19).add(perms_20_24);
  var total_ct = 25;
  var output_scale = perms_0_4.projection().nominalScale();
}

// calculate p-value map
var pVals = allPerms.divide(total_ct).rename('p_val');

// calculate binary mask (where 1s are valid pixels, 0s are masked out)
var sigMask = pVals.lte(maxpVal).rename('signif');
Map.addLayer(sigMask);

// show asnych map, for comparison
var asynch = ee.Image('users/drewhart/NIRv_global_asynch');
var palettes = require('users/gena/packages:palettes');
Map.addLayer(asynch.select('b3'),
             {min:-0.00002, max:0.0002,
               palette: palettes.matplotlib.magma[7], opacity:0.75}, 'asynchrony');
               

// export the significance mask
var roi =
  ee.Geometry.Polygon(
      [[[-179, 70],
        [-179, -70],
        [179, -70],
        [179, 70],
        [-179,70]]], null, false);
Export.image.toAsset({
  image: sigMask.select('signif'),
  description: ee.String(dataset).cat('_permutation_significance_mask').getInfo(),
  assetId: ee.String(dataset).cat('_permutation_significance_mask').getInfo(),
  region: roi,
  scale: output_scale.getInfo(),
  maxPixels: 6e9,
  pyramidingPolicy: {
    '.default': 'mean',
}});

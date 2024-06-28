// load params
var params = require('users/drewhart/seasonality/:params.js');

// combine all permutation fail-count maps
var perms_0_1 = ee.Image('users/drewhart/NIRv_permutations_seeds0-1');
var perms_2_3 = ee.Image('users/drewhart/NIRv_permutations_seeds2-3');
var perms_4_5 = ee.Image('users/drewhart/NIRv_permutations_seeds4-5');
var perms_6_7 = ee.Image('users/drewhart/NIRv_permutations_seeds6-7');
var perms_8_9 = ee.Image('users/drewhart/NIRv_permutations_seeds8-9');
var perms_10_11 = ee.Image('users/drewhart/NIRv_permutations_seeds10-11');
var perms_12_13 = ee.Image('users/drewhart/NIRv_permutations_seeds12-13');
var perms_14_15 = ee.Image('users/drewhart/NIRv_permutations_seeds14-15');
var perms_16_17 = ee.Image('users/drewhart/NIRv_permutations_seeds16-17');
var perms_18_19 = ee.Image('users/drewhart/NIRv_permutations_seeds18-19');
var allPerms = perms_0_1
          .add(perms_2_3)
          .add(perms_4_5)
          .add(perms_6_7)
          .add(perms_8_9)
          .add(perms_10_11)
          .add(perms_12_13)
          .add(perms_14_15)
          .add(perms_16_17)
          .add(perms_18_19);

// calculate p-value map
var total_ct = params.permNIts;
var pVals = allPerms.divide(total_ct).rename('p_val');

// calculate binary mask (where 1s are valid pixels, 0s are masked out)
var sigMask = pVals.lte(params.permMaxPVal).rename('signif');


if (params.map_intermediates){
  Map.addLayer(allPerms, {min: 0, max: total_ct}, 'allPerms');
  Map.addLayer(pVals, {min: 0, max: 1}, 'pVals');
  Map.addLayer(sigMask, {min: 0, max: 1}, 'sigMask');
}


// export the significance mask
Export.image.toAsset({
  image: sigMask.select('signif'),
  description: 'NIRv_permutation_significance_mask',
  assetId: 'NIRv_permutation_significance_mask',
  region: params.roi,
  scale: perms_0_1.projection().nominalScale().getInfo(),
  maxPixels: params.maxPixels,
  pyramidingPolicy: 'mean',
});

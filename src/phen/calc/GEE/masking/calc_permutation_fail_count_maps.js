// calculates a map of the number of times each cell's R2 in a randomly permuted time series
// is greater than that cell's R2 in the unpermuted (i.e., true) time series
// (i.e., the number of times a cell 'fails' the permutation test)
// results will be added together then divided by the total
// number of tests run to develop a global p-value map, for data masking

// NOTE: MUST BE RUN MULTIPLE TIMES, each time updating the permNIts and permSeedStart parameters
//       in params.js. Each time will produce a new GEE asset with permNIts additional permutations
//       captured in it. When enough permutations have been run in total then use
//       calc_permutation_signif_mask.js to summarize them all into a single mask layer.


// load parameters
var params = require('users/drewhart/seasonality/:params.js');


// get requisite functions for data-loading
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
// for significance-testing
var runPerms = require('users/drewhart/seasonality/:fft/run_permutation_tests.js');


// read in the dataset and get its band name, projection, and scale
// NOTE: not doing BRDF albedo QA filtering, to hopefully get around "Execution failed; out of memory" error
var NIRv = gDat.getNIRvData(false,
                            params.maxBRDFAlbQualVal,
                            params.maskWaterSeparately,
                            params.permNIRvNYears,
                            params.permNIRvEndYear,
                            params.permNIRvDayStep,
                            params.NIRvClampToMinPos,
                            params.NIRvUnmaskToMinPos,
                            params.NIRvMaskLteZero
                           );
var dependent = 'NIRv';
var proj = NIRv.first().projection();
var scale = NIRv.first().projection().nominalScale();

// loop over batches of size params.permBatchSize,
// up to and including params.permNIts-1 (because we start from 0)
// as the final seed number,
// and run a batch task for each
for (var seedN = 0; seedN <= params.permNIts - params.permBatchSize; seedN += params.permBatchSize){
  var perms = runPerms.runPermutations(params.harmonics,
                                       params.bandsForViz,
                                       NIRv,
                                       dependent,
                                       params.permBatchSize,
                                       params.permMaxPVal,
                                       seedN
                                      );
  
  // set task and asset names
  var taskAndAssetBaseName = ee.String('NIRv_permutations_seeds');
  var taskAndAssetName = ee.String(taskAndAssetBaseName)
    .cat(ee.String(ee.Number(seedN).toInt()))
    .cat('-')
    .cat(ee.String(ee.Number(seedN + params.permBatchSize - 1).toInt()))
    .getInfo();
  
  
  // export as an asset
  Export.image.toAsset({
    image: perms.select('fail_counts'),
    description: taskAndAssetName,
    assetId: taskAndAssetName,
    region: params.roi,
    scale: scale.getInfo(),
    maxPixels: params.maxPixels,
    pyramidingPolicy: 'mean'
  });

}
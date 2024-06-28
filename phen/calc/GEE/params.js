//---------------------------------------------------------------------------------------------------
// DATA TO PROCESS
//---------------------------------------------------------------------------------------------------
// dataset to use (String that will be later used to set the dataset ImgColl and the dependent var)
exports.datasetName = 'NIRv';
// exports.datasetName = 'SIF';
// exports.datasetName = 'MODISCloud';
// exports.datasetName = 'TerraClimate';

// variable to use, if using TerraClimate data
exports.climateVar = 'pr';
// exports.climateVar = 'def';
// exports.climateVar = 'tmmn';
// exports.climateVar = 'tmmx';



//---------------------------------------------------------------------------------------------------
// BEHAVIORAL
//---------------------------------------------------------------------------------------------------
// params that control what functionalities the scripts should use
exports.verbose = true;

exports.map_intermediates = true;

exports.map = true;

exports.map_result = false;

exports.custom_palette = true;

exports.export_result = true;




//---------------------------------------------------------------------------------------------------
// DATA READING
//---------------------------------------------------------------------------------------------------
exports.maskingMode = 'default';
// exports.maskingMode = 'strict';

// whether or not to mask pixels with low quality BRDF albedo correction
// (only needs to be set to true if maxBRDFAlbQualVal is < 3, as I have confirmed
// that are zero 4s in our entire global time series, and it also seems, from the DAAC specifications
// for both products, as though 4s in MCD43A2 get mapped to 255s in MCD43A4 by MODIS folks,
// and thus get masked out anyhow)
exports.maskLowBRDFAlbQual = false;

// max quality value to allow for the MODIS BRDF and albedo QAQC flags
exports.maxBRDFAlbQualVal = 3;

// whether to mask water (separately from the masking that will be done in LULC filtering)
exports.maskWaterSeparately = false;

// number of years of MODIS NIRv data to use
exports.NIRvNYears = 20;

// last full year of MODIS NIRv data to use
exports.NIRvEndYear = 2021;

// step size between days of MODIS NIRv data to use
exports.NIRvDayStep = 4;

// whether to clamp NIRv values <=0 to the minimum positive observed value
// (to backfill and prevent data dropout in some parts of boreal regions)
exports.NIRvClampToMinPos = true;

// whether to unmask missing vals to the minimum positive NIRv value observed at those pixels
// (to explore places with extended seasonal gaps; not used in final analysis)
exports.NIRvUnmaskToMinPos = false;

// whether to mask out values <= zero
// NOTE: will always be true when NIRv data itself is actually being used;
//       only false when calculating monthly proportional data availability, to cut mem usage
exports.NIRvMaskLteZero = true;




//---------------------------------------------------------------------------------------------------
// PERMUTATION TESTING
//---------------------------------------------------------------------------------------------------
// total number of permutations
// NOTE: SHOULD NEVER EXCEED 232, OTHERWISE UNSIGNED 8-BIT INTEGER WILL FAIL!
exports.permNIts = 20;

// how many permutations to do per batch
exports.permBatchSize = 2;

// max P-value beyond which to mask out a pixel (i.e., alpha value)
// NOTE: effective alpha will depend on the smallest value possible for the
//       the number of permutation tests being used (i.e., on the value
//       of permNIts), which will be 0.05 if we enforce a rule
//       that none of our 20 iterations can have 'failed'
exports.permMaxPVal = 0.001;

// because of computational constraints on GEE
// (definitely not designed for this sort of analysis),
// use a longer day step and shorter total time series for the NIRv data
// than what is defined in the DATA READING AND MASKING section above
// (which, NOTE, will only make the permutation-test results more conservative)
exports.permNIRvDayStep = 7;
exports.permNIRvNYears = 10;
exports.permNIRvEndYear = 2016;




//---------------------------------------------------------------------------------------------------
// OTHER MASKING
//---------------------------------------------------------------------------------------------------
// minimum mean monthly proportion of valid days
// (below which to mask out a pixel for not having adequate monthly coverage)
exports.minMonthlyPropDataAvail = 0.1;

// minimum Pielou evenness of MODIS time series 
exports.minTSEvenness = 0.8;

// min percent of non-null values needed in the time series 
// (less than this will be filtered)
exports.minPctTSDataAvailability = 0.5;

// whether to plot and export individual masks
exports.plotAndExportMasks = true;




//---------------------------------------------------------------------------------------------------
// LSP MODELING
//---------------------------------------------------------------------------------------------------
// number of harmonics to use in harmonic regression (2 means one annual, one semiannual)
var nFreq = 2;
exports.harmonics = ee.List.sequence(1, nFreq); 

// bands to use for calculating phase and amplitude
exports.bandsForViz = ee.List(['sin_1', 'cos_1']);






//---------------------------------------------------------------------------------------------------
// FILE EXPORT
//---------------------------------------------------------------------------------------------------

// DETH: changing from 256x256 patches to 300x300
//       because that evenly breaks up the [(-165,60), (180,-60)] quasi-global study area
//       into 184 even-sized patches (23 cols x 8 rows, 
//       with each patch being 476*476 after adding the 176-cell margins);
//       each patch should be approx. 4.3 MB, so if I set the maxFileSize to 5_000_000 bytes then
//       I should get 184 4.3MB files...

// global ROI
exports.roi =
  ee.Geometry.Polygon(
      [[[-165, 75],
        [-165, -60],
        [180, -60],
        [180, 75],
        [-165, 75]]], null, false);

// whether to run chunked, separate export tasks
exports.in_chunks = false;

// max file size in bytes
exports.maxFileSize = 10000000;

// file format can be 'TFRecord', 'Asset', or 'GeoTIFF'
exports.fileFormat = 'TFRecord';

//square-patch size, in pixels
exports.patchSize = 300;

// maxPixels argument fed to export.Image.toDrive()
exports.maxPixels = 4000000000;

// shard size to use for chunked processing (note GeoTIFF file dimensions must be multiples of this!)
exports.shardSize = 256;

// GeoTIFF file dimensions, used used for all mask map files
exports.fileDimensions = [2816, 6912];

// max neighbor distance to be used in calculating asynchrony (in meters)
// (NOTE: actual value would be 150000 m, so this is a comfortable overestimate)
exports.maxNeighDist = ee.Number(180000);

// factor by which to multiply the calculated min-needed kernel size,
// to ensure adequate kernel size (in pixels) to really get all neighbors within max neigh dist (in m)
exports.kernelSizeBufferFactor = 1.8;

// default value to use to fill missing data
exports.defaultValue = -9999;

// abbreviated name of region enclosed by the bounding box (within which outputs will be generated)
exports.regionAbbrev = ee.String('global');

// Google Drive folder where all main and mask file exports should be saved 
exports.mainFolder = 'LSP_outputs_from_GEE';
exports.maskFolder = 'LSP_mask_outputs_from_GEE';




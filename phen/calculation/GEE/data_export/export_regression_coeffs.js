//-------
// PARAMS
//-------


//-------
// main output controls:

// dataset to use (String that will be later used to set the dataset ImgColl and the dependent var)
//var datasetName = 'SIF';
var datasetName = 'NIRvP';

var map_result = false;
var export_result = true;
var in_chunks = false;
var maxFileSize = 7000000; // in bytes
// DETH: changing from 256x256 patches to 300x300
//       because that evenly breaks up the [(-165,60), (180,-60)] quasi-global study area
//       into 184 even-sized patches (23 cols x 8 rows, 
//       with each patch being 476*476 after adding the 176-cell margins);
//       each patch should be approx. 4.3 MB, so if I set the maxFileSize to 5_000_000 bytes then
//       I should get 184 4.3MB files...
var patchSize = 300; //square-patch size, in pixels
var maxPixels = 20000000; // maxPixels argument fed to export.Image.toDrive()
var fileFormat = 'TFRecord'; // 'GeoTIFF'
// max neighbor distance to be used in calculating asynchrony (in meters)
var maxNeighDist = ee.Number(180000);
// factor by which to multiply the calculated min-needed kernel size,
// to ensure adequate kernel size (in pixels) to really get all neighbors within max neigh dist (in m)
var kernelSizeBufferFactor = 1.8;
// abbreviated name of region enclosed by the bounding box (within which outputs will be generated)
var regionAbbrev = ee.String('global');
var regionFolder = ee.String('asynch_outputs_from_GEE');
var globalFolder = ee.String('asynch_outputs_from_GEE');

//-------

// number of harmonics to use in harmonic regression (2 means one annual, one semiannual)
var nFreq = 2;
var harmonics = ee.List.sequence(1, nFreq); 
// and bands to use for calculating phase and amplitude
var bandsForViz = ee.List(['sin_1', 'cos_1']);

// min percent of observ (over time) at a pix to cause it to be masked
var minWaterOccur = 80;

// min percent valid land cover among the land-cover pixels that fall inside
// each input dataset (NIRv, SIF) pixel (pixels with less than this will be filtered out)
var minPctValidLC = 0.75;

// min perecent of non-null values needed in the time series 
// (less than this will be filtered)
var minPctTSLen = 0.5;


// NOTE: COMMENT OUT PERMUTATION TEST FOR NOW
/*
// number of permutations,
// max p-value to use for the significance test,
// and number to use to start seq of integers that serve as seeds for each iteration
var nIts = 10;
var maxPVal = 0.05;
var seedStart = 0;
*/

// min percent of valid neighbors that a pixel must have
// in order to be included in the output asynchrony map
var minPctValidNeigh = 0.5;


//---------
// REQUIRES
//---------

// data-loading
var gDat = require('users/drewhart/seasonality:io/get_data.js');
// filtering
var maskLC = require('users/drewhart/seasonality:masking/mask_landcover.js');
var maskShortTS = require('users/drewhart/seasonality:masking/mask_short_time_series.js');
var getSigMask = require('users/drewhart/seasonality:masking/get_significance_mask.js');
// harmonic regression
var hReg = require('users/drewhart/seasonality:fft/harmonic_regression.js');
// annual fitted time series
var fitYr = require('users/drewhart/seasonality:fft/fit_single_year_daily_values.js');
// asynchrony
var asynch = require('users/drewhart/seasonality:asynchrony/calc_asynchrony.js');
// map stuff
var cbar = require('users/drewhart/seasonality:viz/make_colorbar.js');
var palettes = require('users/gena/packages:palettes'); 



//--------
// IMPORTS
//--------

// read in the dataset needed
if (datasetName === 'NIRvP'){
  // NIRvP, its dependent variable band name, proj, scale, and nNeigh
  var dataset = gDat.getNIRvPData(minWaterOccur);
  var dependent = 'NIRvP';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
  //var nNeigh = maxNeighDist.divide(scale).ceil(); // num neighs for asynch calc
  
  // or SIF, its dependent variable band name, proj, scale, and nNeigh
} else if (datasetName === 'SIF'){
  var dataset = gDat.getSIFData(minWaterOccur);
  var dependent = 'b1';
  var proj = dataset.first().projection();
  var scale = dataset.first().projection().nominalScale();
  //var nNeigh = maxNeighDist.divide(scale).ceil(); // num neighs for asynch calc
}




//==================
// MASKING/FILTERING
//==================

//------------------
// filter land cover
//------------------
// NOTE: this will also reduce the temporal extent of the data
//       if the data's original temporal extent falls outside that of the LC data!
var LCMasked = maskLC.maskLandCover(dataset, minPctValidLC, proj);


//-------------------
// filter missingness
//-------------------
var LCTSMasked = maskShortTS.maskShortTimeSeries(LCMasked, dependent, minPctTSLen);


// calculate regression with nFreq freqs
var reg = ee.ImageCollection(hReg.calcHarmonicRegression(LCTSMasked, dependent,
                                                         harmonics, bandsForViz));
                                                         


// NOTE: COMMENT OUT PERMUTATION TESTS FOR NOW:
/*
//---------------------------------------------------
// filter by significance of permutation test of R^2s
//---------------------------------------------------

// run permutation test and get mask
// NOTE: leaving permutation test function as returning a mask, not a masked ImgColl,
//       because I may wind up wanting to produce and export the mask (or multiples of
//       them) separately, then read them in here to use them
var sigMask = getSigMask.runPermutationTest(LCTSMasked, dependent, nIts, maxPVal, seedStart);
var LCTSSigMasked = LCTSMasked.map(function(img){
  return img.updateMask(sigMask)});
*/


if (map_result){
  //-------
  // map it
  //-------
  Map.setOptions('TERRAIN');
  Map.addLayer(reg.first().select('phase'),
               {palette: palettes.kovesi.cyclic_mygbm_30_95_c78_s25[7], 
                min:0, max:1, opacity: 0.8}, 'phase');
  Map.add(cbar.makeLegend('phase', 'Jan 1', 'Jul 1', 'Dec 31',
          palettes.kovesi.cyclic_mygbm_30_95_c78_s25[7], 'bottom-center'));
}


// export the regression coeffs
if (export_result){
  if (in_chunks){
    //----------
    // EXPORT IT
    //----------
    
    var lon_chunks = [[-170, -135], [-140, -105], [-110, -75], [-80, -45],
                      [-50, -15], [-20, 15], [10, 45], [40, 75], [70, 105],
                      [100, 135], [130, 165], [160, 179]];
    var lat_chunks = [[60, 20], [25, -25], [-20, -60]];
    for (var lons_i=0; lons_i < lon_chunks.length; lons_i++) {
      for (var lats_i=0; lats_i < lat_chunks.length; lats_i++) {
        var lon1 = lon_chunks[lons_i][0];
        var lon2 = lon_chunks[lons_i][1];
        var lat1 = lat_chunks[lats_i][0];
        var lat2 = lat_chunks[lats_i][1];
        var roi = ee.Geometry.Polygon([[[lon1, lat1],
                                        [lon1, lat2],
                                        [lon2, lat2],
                                        [lon2, lat1],
                                        [lon1, lat1]]], null, false);
    
    var output = reg
      .first()
      .select(['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']);
    
    var scale = dataset
      .first()
      .projection()
      .nominalScale()
      .getInfo();
    
    // calculate the needed kernelSize in order to get adequate patch-overlap margins that allow for
    // calculation of asynchrony all the way out to the neighborhood-radius distance
    // NOTE: multiply by 2 because overlap margin is of size kernelSize/2
    // DETH: 10-09-21: JUST REALIZED! this is not providing adequate marginal overlap
    //                 at higher latitudes because I'm using the nominal scale in meters
    //                 **AT THE POINT OF TRUE SCALE**, which is near the equator, not toward the poles!
    //                 changing this instead to hard-code for the minimum per-cell distance in our dataset
    //                 (i.e. the longitudinal distance across a cell all the way up/down at 60 N/S,
    //                 given that we're working in a fixed res of ~0.05 degrees),
    //                 which according to https://en.wikipedia.org/wiki/Longitude
    //                 is ~55.8 km per degree * 0.05 degrees = ~2.79km per cell,
    //                 then rounding down further to 2.25km, to ensure that we take more than adequate margin sizes...
    var kernelSize = maxNeighDist
      //.divide(scale)
      .divide(ee.Number(2.25*1000))
      .multiply(kernelSizeBufferFactor)
      .multiply(2)
      .ceil()
      .getInfo();
    // add 1 to the kernelSize if it's odd (since it needs to be evenly divisible by 2)
    if (kernelSize % 2 == 1){
      var kernelSize = ee.Number(kernelSize).add(1).getInfo();
    }
    print('kernelSize', kernelSize, 'scale', scale, 'maxNeighDist', maxNeighDist);
     
   
    var lons_lats_str =  ee.String('_')
      .cat(ee.Number(lon1).toInt().format('%03d'))
      .cat('_')
      .cat(ee.Number(lat1).toInt().format('%02d'))
      .cat('_to_')
      .cat(ee.Number(lon2).toInt().format('%03d'))
      .cat('_')
      .cat(ee.Number(lat2).toInt().format('%02d'));
    var taskAndAssetName = regionAbbrev.cat('_coeffs_').cat(datasetName).cat(lons_lats_str).getInfo();
  
    // Export the image, specifying scale and region.
    if (fileFormat == 'TFRecord'){
    Export.image.toDrive({
      image: output,
      description: taskAndAssetName,
      //folder: regionFolder.getInfo(),
      folder: globalFolder.getInfo(),
      scale: dataset.first().projection().nominalScale().getInfo(),
      region: roi,
      // use the file format provided at the top of the script
      fileFormat: fileFormat,
      maxPixels: maxPixels,
      formatOptions: {
        maxFileSize: maxFileSize,
        defaultValue: -9999,
        // NOTE: It appears that both patchDimensions AND kernelSize need to be expressed in
        // pixel count, not in meters...
        patchDimensions: [patchSize, patchSize],
        // NOTE: kernelSize should be such that kernelSize/2 > maxDist used in asynch calc;
        //       because in pixel count, not meters, if I need >300km overlap around edge
        //       (because that's my target maxDist that defines the neighborhood within which
        //       asynchrony is calculated), that's >300,000m = 1/2 kernel Size,
        //       so 600,000m = kernel size in meters,
        //       so 600,000/6000 (nominal scale in meters, rounded up from 5565.974539664152m)
        //       = 100 pixels, so that's my kernelSize that will ensure an adequate margin to
        //       calculate asynchrony from each patch without crosstalk between patches
        kernelSize: [kernelSize, kernelSize]
       }
     });
    }
    else if (fileFormat == 'GeoTIFF'){
    Export.image.toDrive({
      image: output,
      description: taskAndAssetName,
      //folder: regionFolder.getInfo(),
      folder: globalFolder.getInfo(),
      scale: dataset.first().projection().nominalScale().getInfo(),
      region: roi,
      // use the file format provided at the top of the script
      fileFormat: fileFormat,
      maxPixels: maxPixels,
      formatOptions: {
        defaultValue: -9999,
        // NOTE: It appears that both patchDimensions AND kernelSize need to be expressed in
        // pixel count, not in meters...
        patchDimensions: [patchSize, patchSize],
        // NOTE: kernelSize should be such that kernelSize/2 > maxDist used in asynch calc;
        //       because in pixel count, not meters, if I need >300km overlap around edge
        //       (because that's my target maxDist that defines the neighborhood within which
        //       asynchrony is calculated), that's >300,000m = 1/2 kernel Size,
        //       so 600,000m = kernel size in meters,
        //       so 600,000/6000 (nominal scale in meters, rounded up from 5565.974539664152m)
        //       = 100 pixels, so that's my kernelSize that will ensure an adequate margin to
        //       calculate asynchrony from each patch without crosstalk between patches
        kernelSize: [kernelSize, kernelSize]
       }
     });
     }
    }
   }
  }
else {
  //----------
  // EXPORT IT
  //----------
  var roi =
    ee.Geometry.Polygon(
        [[[-165, 60],
          [-165, -60],
          [180, -60],
          [180, 60],
          [-165,60]]], null, false);
    
    var output = reg
      .first()
      .select(['constant', 'sin_1', 'cos_1', 'sin_2', 'cos_2']);
    
    var scale = dataset
      .first()
      .projection()
      .nominalScale()
      .getInfo();
    
    // calculate the needed kernelSize in order to get adequate patch-overlap margins that allow for
    // calculation of asynchrony all the way out to the neighborhood-radius distance
    // NOTE: multiply by 2 because overlap margin is of size kernelSize/2
    // DETH: 10-09-21: JUST REALIZED! this is not providing adequate marginal overlap
    //                 at higher latitudes because I'm using the nominal scale in meters
    //                 **AT THE POINT OF TRUE SCALE**, which is near the equator, not toward the poles!
    //                 changing this instead to hard-code for the minimum per-cell distance in our dataset
    //                 (i.e. the longitudinal distance across a cell all the way up/down at 60 N/S,
    //                 given that we're working in a fixed res of ~0.05 degrees),
    //                 which according to https://en.wikipedia.org/wiki/Longitude
    //                 is ~55.8 km per degree * 0.05 degrees = ~2.79km per cell,
    //                 then rounding down further to 2.25km, to ensure that we take more than adequate margin sizes...
    var kernelSize = maxNeighDist
      //.divide(scale)
      .divide(ee.Number(2.25*1000))
      .multiply(kernelSizeBufferFactor)
      .multiply(2)
      .ceil()
      .getInfo();
    // add 1 to the kernelSize if it's odd (since it needs to be evenly divisible by 2)
    if (kernelSize % 2 == 1){
      var kernelSize = ee.Number(kernelSize).add(1).getInfo();
    }
    print('kernelSize', kernelSize, 'scale', scale, 'maxNeighDist', maxNeighDist);
    
    var taskAndAssetName = regionAbbrev.cat('_coeffs_').cat(datasetName).getInfo();
  
    // Export the image, specifying scale and region.
    if (fileFormat == 'TFRecord'){
    Export.image.toDrive({
      image: output,
      description: taskAndAssetName,
      //folder: regionFolder.getInfo(),
      folder: globalFolder.getInfo(),
      scale: dataset.first().projection().nominalScale().getInfo(),
      region: roi,
      // use the file format provided at the top of the script
      fileFormat: fileFormat,
      maxPixels: maxPixels,
      formatOptions: {
        maxFileSize: maxFileSize,
        defaultValue: -9999,
        // NOTE: It appears that both patchDimensions AND kernelSize need to be expressed in
        // pixel count, not in meters...
        patchDimensions: [patchSize, patchSize],
        // NOTE: kernelSize should be such that kernelSize/2 > maxDist used in asynch calc;
        //       because in pixel count, not meters, if I need >300km overlap around edge
        //       (because that's my target maxDist that defines the neighborhood within which
        //       asynchrony is calculated), that's >300,000m = 1/2 kernel Size,
        //       so 600,000m = kernel size in meters,
        //       so 600,000/6000 (nominal scale in meters, rounded up from 5565.974539664152m)
        //       = 100 pixels, so that's my kernelSize that will ensure an adequate margin to
        //       calculate asynchrony from each patch without crosstalk between patches
        kernelSize: [kernelSize, kernelSize]
       }
     });
    }
    else if (fileFormat == 'GeoTIFF'){
    Export.image.toDrive({
      image: output,
      description: taskAndAssetName,
      //folder: regionFolder.getInfo(),
      folder: globalFolder.getInfo(),
      scale: dataset.first().projection().nominalScale().getInfo(),
      region: roi,
      // use the file format provided at the top of the script
      fileFormat: fileFormat,
      maxPixels: maxPixels,
      formatOptions: {
        defaultValue: -9999,
        // NOTE: It appears that both patchDimensions AND kernelSize need to be expressed in
        // pixel count, not in meters...
        patchDimensions: [patchSize, patchSize],
        // NOTE: kernelSize should be such that kernelSize/2 > maxDist used in asynch calc;
        //       because in pixel count, not meters, if I need >300km overlap around edge
        //       (because that's my target maxDist that defines the neighborhood within which
        //       asynchrony is calculated), that's >300,000m = 1/2 kernel Size,
        //       so 600,000m = kernel size in meters,
        //       so 600,000/6000 (nominal scale in meters, rounded up from 5565.974539664152m)
        //       = 100 pixels, so that's my kernelSize that will ensure an adequate margin to
        //       calculate asynchrony from each patch without crosstalk between patches
        kernelSize: [kernelSize, kernelSize]
       }
     });
    }
  }
}
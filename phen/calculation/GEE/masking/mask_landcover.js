/* Final form of the landcover masking code,
   combined from Sarah's and Drew's work

  EXPORTS:
    maskLandCover(imgColl, minPctValidLC, proj)
      imgColl: the ImageCollection for which to produce the mask
      minPctValidLC: the min percent valid LC pixels within
                     an input ImageCollection pixel that is 
                     required for that pixel to be included 
                     in the output masked ImageCollection
      proj: the projection to which to reproject the LC before masking with it
      
      returns: the masked ImageCollection
      
      NOTE: the returned ImageCollection can be shorter in time extent
            because dates outside of the date range of the LC data are being dropped
*/
   

// load the MODIS landcover dataset (LC classification Type 1)
var MODISLC = ee.ImageCollection("MODIS/006/MCD12Q1")
  // NOTE: using LC_Type1 rather than Type2 because it's an international standard
  // classification scheme, as explained by the data's User Guide
  .select('LC_Type1')
  // and add 'year' as a property, to use when later masking images by their year
  .map(function(img){
    return img.set('year', ee.Date(ee.Image(img).date()).get('year'))});

// get the date range covered by the LC data, so that the SIF and NIRv datasets can
// be filtered to it
// (NOTE: at time of writing, 08/14/20, this only matters for the NIRv)
var firstLCYear = MODISLC.first().date().get('year');
var lastLCYear = ee.Image(MODISLC.toList(1, MODISLC.size().subtract(1)).get(0)).date().get('year');

// function for masking invalid landcover out of an input ImageCollection
exports.maskLandCover = function(maskingMode, imgColl, minPctValidLC, proj, preLCMask) {
  // filter the imgColl to only the dates covered by our LC data
  print('dates before', 
        imgColl.first().date(), 
        ee.Image(imgColl.toList(1, imgColl.size().subtract(1)).get(0)).date());
  // get max LC class number, based on masking mode
  if (maskingMode == 'strict'){
    var maxLCClass = 12;
  } else if (maskingMode == 'default'){
    var maxLCClass = 15;
  }
  var imgCollFilteredDates = imgColl.filter(ee.Filter.calendarRange(firstLCYear, lastLCYear, 'year'));
  print('dates after',
        imgCollFilteredDates.first().date(),
        ee.Image(imgCollFilteredDates.toList(1, imgCollFilteredDates.size().subtract(1)).get(0)).date());
  var imgCollFilteredLC = imgCollFilteredDates.map(function(img){
    // get the current image's year
    var year = img.date().get('year');
    // then grab the MODIS LC image that matches that year
    var LCImg = MODISLC.filter(ee.Filter.eq('year', year)).first();
    // then use that LC image to tag for masking LC classes >= 12 (for strict mask)
    // or >= 15 (for default). These include: 
    //      12: cropland
    //      13: urban
    //      14: crop-natural mosaic
    //      15: permanent snow/ice
    //      16: barren
    //      17: permanent water bodies
    var LCMask =　LCImg.lt(maxLCClass)
      // reproject to the CRS of the input ImageCollection, 
      // using a mean to reduceResolution if necessary
      // (so that the resulting values represent the percent of valid LC pixels contained
      //  within an output ImageCollection pixel, i.e. the percent in LC classes <12)
      .reduceResolution({
        reducer: ee.Reducer.mean(),
        maxPixels: 2000})
      .reproject({crs: proj})
      // create a mask in which only cells with >= minPctValidLC remain unmasked
      .gte(minPctValidLC);
    // mask the image and return it
    return img.updateMask(LCMask.and(preLCMask));
  });
  return imgCollFilteredLC;
};

      



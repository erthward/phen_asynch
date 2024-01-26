//////////////////////
//LAUREN'S FUNCTION//
////////////////////

exports.maskShortTimeSeries = function(imgColl, bandName, minPctTSLen){
  // function to reduce an image collection using count that filters out all null values and 
  /* values less than a set minimum time-series length ('minTSLength') 
      KEY:
        imgColl     : takes an image collection
        bandName    : the name of the band you want the function to operate on;
                      must be formatted as such with '' on either side: 'imgColl_count'; 
                      e.g. imgColl = SIF, bandName = 'b1'
        minPctTSLen : the desired minumum time-series length, as a percent of all TS steps
                      that must have null values; takes in a fraction between 0 and 1
  */
  // set the min and max time-series lengths based on the minPctTSLen
  var maxTSLen = ee.Number(imgColl.size());
  var minTSLen = ee.Number(maxTSLen.multiply(minPctTSLen)).toInt();
  // create an image that has 1's at all pixels with count > minTSLength
  var minTSLenMask = imgColl
    .select(bandName)
    .reduce(ee.Reducer.count())
    .gte(minTSLen);
  // mask the ImageCollection with that mask
  //var maskedImgColl = imgColl.map(function(img){
  //  return img.updateMask(minTSLenMask)});
  //return maskedImgColl; 
  return minTSLenMask;
};
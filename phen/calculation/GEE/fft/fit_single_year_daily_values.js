/* EXPORTS:

fitSingleYearDailyValues = function(harmRegResults, proj)
    harmRegResults: an Image containing the summary results of an annual-and-semiannual
                    harmonic regression model (the first Image in the ImageCollection
                    output by our harmonic regression function)
    proj: the projection in which to calculate and return the year of daily fitted values
  
*/



/* TODO:

- how to correctly assign system:time_start to each image, rather than just num_index???

- why do some sites wind up with values in the negatives???

*/



//-------
// params
//-------


var pi = ee.Number(3.141592653589793);
var msecPerDay = ee.Number(86400000);

/*
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
var SIF = gDat.getSIFData(50);
var SIFProj = SIF.first().projection();

// calculate regression with n_freq freqs
var hReg = require('users/drewhart/seasonality/:fft/harmonic_regression.js');
var n_freq = 2;
var harmonics = ee.List.sequence(1, n_freq); 
var dependent = 'b1';
var bandsForViz = ee.List(['sin_1', 'cos_1']);
var reg = ee.ImageCollection(hReg.calcHarmonicRegression(SIF, dependent,
                                                         harmonics, bandsForViz));
*/



//----------
// functions
//----------

// function to create an ImageCollection containing 365 Images,
// one for each day of the year, with each one having bands for
// the day expressed as the sin or cos of its radians on either
// an annual or semiannual frequency
// NOTE: creating a set of fitted values for these dates will not only give me
//       an equal-length set of values for the NIRv and SIF datasets from which
//       to calculate asynchrony, but will also detrend those datasets by dropping
//       the linear-time covariate
var makeSingleYearHarmonicImgColl = function(datasetProj){
  // get the integer days, 0 to 364
  var days = ee.List.sequence({start:0, end:364, step:1});
  // re-express those as radians in both the annual and semiannual freqs
  var annualRads = days.map(function(n){
    return ee.Number(n).divide(365).multiply(pi.multiply(2))});
  var semianRads = days.map(function(n){
    return ee.Number(n).divide(365).multiply(pi.multiply(4))});
  // create ImageCollections of constant images of the sines and cosines
  // of those radian-transformed values
  var annualSin = ee.ImageCollection(annualRads.map(function(n){
    return ee.Image.constant(n).sin().rename('ann_sin')}));
  var annualCos = ee.ImageCollection(annualRads.map(function(n){
    return ee.Image.constant(n).cos().rename('ann_cos')}));
  var semianSin = ee.ImageCollection(semianRads.map(function(n){
    return ee.Image.constant(n).sin().rename('sem_sin')}));
  var semianCos = ee.ImageCollection(semianRads.map(function(n){
    return ee.Image.constant(n).cos().rename('sem_cos')}));
  // combine all into one ImageCollection
  var year = ee.ImageCollection(ee.ImageCollection(ee.ImageCollection(annualSin
    .combine(annualCos))
    .combine(semianSin))
    .combine(semianCos))
  // reproject to the datasetProj
    .map(function(img){return img.reproject(datasetProj)});
  var yearWithIndex = year.map(function(img){
    return img.set({'num_index': ee.Number.parse(ee.Image(img).get('system:index')),
      // 'system:time_start': ee.Date(ee.Number.parse(ee.Image(img).get('system:index')).multiply(msecPerDay))})})
    })})
    .sort('num_index');
  return yearWithIndex;
};


// function that takes the first image from a harmonic regression ImgColl
// (the image which contains all the coefficients and so on)
// and a year's worth of annual and semiannual-freq sine and cosine harmonics,
// and returns a new 365-length ImgColl of the fitted values for each day of the year
exports.fitSingleYearDailyValues = function(harmRegResults, proj){
  var year = makeSingleYearHarmonicImgColl(proj);
  var fittedYear = year.map(function(img){
    return harmRegResults.select('constant')
      .add(harmRegResults.select('sin_1').multiply(img.select('ann_sin')))
      .add(harmRegResults.select('cos_1').multiply(img.select('ann_cos')))
      .add(harmRegResults.select('sin_2').multiply(img.select('sem_sin')))
      .add(harmRegResults.select('cos_2').multiply(img.select('sem_cos')))
      .rename('fitted')
      // copy the properties, to carry over the num_index property I created
      .copyProperties(img)});
  return fittedYear;
};


 
//var year = makeSingleYearHarmonicImgColl(SIFProj);
//print('year', year);
//Map.addLayer(year, {}, 'year');

//var fittedYear = exports.fitSingleYear(reg.first(), SIFProj);
//print('fittedYear', fittedYear);
//Map.addLayer(fittedYear, {}, 'fittedYear');
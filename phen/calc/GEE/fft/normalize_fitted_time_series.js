
//Normalize all pixels in that dataset output by 'fit_single_year_daily_values'
//so that their fitted values fluctuate between 0 and 1 at their respective 
//mins and maxes, returning an ImageCollection

//////////////////////////////////////////////////
//Read in Data, Imports, Paramters, & Functions//
////////////////////////////////////////////////

// read in the SIF data (with 50% min water occurrence)
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
var SIF = gDat.getSIFData(50);
//Read in parameters for 'fitSingleYearDailyValues' function
var SIFProj = SIF.first().projection();

// Read in functions 
var fittedVals = require('users/drewhart/seasonality/:fft/fit_single_year_daily_values.js');
var hReg = require('users/drewhart/seasonality/:fft/harmonic_regression.js');

//Read in parameters of 'calcHarmonicRegression' function
var n_freq = 2;
var harmonics = ee.List.sequence(1, n_freq); 
var dependent = 'b1';
var bandsForViz = ee.List(['sin_1', 'cos_1']);

//
//Make Image Collections from Imported functions
//

// make image collection of regression outputs
var reg = ee.ImageCollection(hReg.calcHarmonicRegression(SIF, dependent,
                                                         harmonics, bandsForViz));
//print("reg", reg);

// get the fitted single year daily values from the haromic regression
// from 'calcHarmonicRegression' function
var fittedSingleYearDVImgColl = ee.ImageCollection(fittedVals.fitSingleYearDailyValues(
                                                                              reg.first(), 
                                                                              SIFProj));
print("fittedSingleYearDVImgColl", fittedSingleYearDVImgColl);

////////////////////////////
//Normalize Fitted Values//
//////////////////////////

//Normalize all pixels in that dataset output of fitted single year 
//daily regression values (fluctuate between 0 and 1 at their respective 
// mins and maxes); returns an ImageCollection
  //KEY: 
      //ImgColl: an ImageCollection of fitted single year daily 
          //regression values from the 'calcHarmonicRegression' function
exports.normalizeFittedVals = function(ImgColl){
  //get the minimum value Image in the ImageCollection
    var min = ImgColl.reduce(ee.Reducer.min());
  //get the maximum value Image in the ImageCollection
    var max = ImgColl.reduce(ee.Reducer.max());
  //compute the numerator for the normalization expression
    var numerator = ImgColl.map(function(image){
      return ee.Image(image.subtract(min))});
  //compute the denominator for the normalization expression
    var denominator = max.subtract(min);
  //map normalization expression across entire ImageCollection
    var normImgColl = numerator.map(function(image){
      return ee.Image(image.divide(denominator))});
  return normImgColl;
};

// TEST THE FUNCTION
//var normTS = NormFittedVals(fittedSingleYearDVImgColl);
//print("normTS", normTS);

// Plot a map of the time series
//Map.addLayer(normTS, 
//    {palette: ['blue'], min: 0, max: 365},
//   'TS');



/////////////////////////
//TEST CODE FOR FUNCTION

//get the minimum value in the image
//var min = fittedSingleYearDVImgColl.reduce(ee.Reducer.min());
//print("min", min);

//get the maximum value in the image 
//var max = fittedSingleYearDVImgColl.reduce(ee.Reducer.max());
//print("max", max);

//compute the numerator for the normalization expression
//var numerator = fittedSingleYearDVImgColl.map(function(image){
//  return ee.Image(image.subtract(min))});
//print("numerator", numerator);
  
//var denominator = max.subtract(min);
//print("denominator", denominator);

//var normImgColl = numerator.map(function(image){
//  return ee.Image(image.divide(denominator))});
//print("normImgColl", normImgColl);


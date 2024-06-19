/*
* Functions for calculating harmonic linear regression.
*  
* STARTING WITH CODE TAKEN FROM: 
* https://gis.stackexchange.com/questions/271454/harmonic-modis-
* trend-in-google-earth-engine
* which I think, in turn, came from Nick Clinton's tutorials in the EDU tab of GEE docs
*/


/* EXPORTS:

calcHarmonicRegression(imgColl, dependent, harmonics, bandsForViz)

  imgColl: an ImageCollection to be used for the harmonic regression
  dependent: the band name to be used as the dependent variable
  harmonics: a List of the harmonic frequencies to be used in the regression 
             (expressed in cycles/year, e.g. 1 = annual, 2 = twice annually...)
  bandsForViz: which sin, cos pair of output bandnames from the regression to be used to calculate the phase
               metric to be visualized, named as ['sin_N', 'cos_N'], where 'N' is the harmonic
               frequency to be used)
  
*/




//////////////////////
// SECONDARY FUNCTIONS
//////////////////////

// Function to get a sequence of band names for harmonic terms.
var constructBandNames = function(base, list) {
  return ee.List(list).map(function(i) {
    return ee.String(base).cat(ee.Number(i).int());
  });
};

// Function to add a time band.
var addDependents = function(img) {
  // Compute time in fractional years since the epoch.
  var unreformattedTime = ee.Date(ee.Number(ee.Image(img).get('t')).multiply(100000));
  var years = unreformattedTime.difference('2000-01-01', 'year');
  var timeRadians = ee.Image(years.multiply(2 * Math.PI)).rename('t');
  var constant = ee.Image(1);
  return ee.Image(img).addBands(constant).addBands(timeRadians.float());
};

// Function to compute the specified number of harmonics
// and add them as bands.  Assumes the time band is present.
var addHarmonics = function(freqs, sinNames, cosNames) {
  return function(img) {
    // Make an image of frequencies.
    var frequencies = ee.Image.constant(freqs);
    // This band should represent time in radians.
    var time = ee.Image(img).select('t');
    // Get the cosine terms.
    var cosines = time.multiply(frequencies).cos().rename(cosNames);
    // Get the sin terms.
    var sines = time.multiply(frequencies).sin().rename(sinNames);
    return ee.Image(img).addBands(cosines).addBands(sines);
  };
};

// function to compute R2s and return as an Image
// (assumes an ImageCollection with Images each having a band of fitted values 
// with the name 'fitted' and a band of original dependent-variable values with
// the name being whatever value the argument `dependent` holds)
var calcR2 = function(imgColl, dependent){
  return imgColl.select('fitted').reduce(ee.Reducer.variance())
          .divide(imgColl.select(dependent).reduce(ee.Reducer.variance()))
          .rename('R2');
};


////////////////
// MAIN FUNCTION
////////////////

// Function for calculating a harmonic regression (fitted-value bands
// for all Images in the ImageCollection)
// and its summary stats (coeffs, resids, neighborhood SD, etc)
// ARGS:
//      imgColl:      the ImageCollection to calculate the harmonic regression for
//      dependent:    the band name of the dependent variable
//      harmonics:    a list of the harmonic freqeuencies to use (1=annual, 2=biannual, etc.)
//      bandsForViz:  a list of the bands to be used to calculate phase, amplitude, for viz
//
exports.calcHarmonicRegression = function(imgColl, dependent, harmonics, bandsForViz){
 
  //------------------------------------
  // Create the harmonic ImageCollection
  //------------------------------------
  
  // Construct lists of names for the harmonic terms.
  var cosNames = constructBandNames('cos_', harmonics);
  var sinNames = constructBandNames('sin_', harmonics);
 // Independent variables.
  var independents = ee.List(['constant', 't'])
    .cat(cosNames).cat(sinNames);
    
  // Add variables
  var harmonicImgColl = ee.ImageCollection(imgColl
    // add the dependent and harmonic variables
    .map(addDependents)
    .map(addHarmonics(harmonics, sinNames, cosNames))
  // and add independents as a property
    .set({'independents': independents}));

  //-------------------------
  // Calculate the regression
  //-------------------------
  
  // Calculate regression reduction and get output (a 4x1 array image)
  var harmonicTrend = harmonicImgColl
    .select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
  //print('harmonicTrend', harmonicTrend);
  
  // Turn the array image into a multi-band image of coefficients.
  var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);
  //print('harmonicTrendCoefficients', harmonicTrendCoefficients);
  
  // calculate the residuals
  var harmonicTrendResiduals = harmonicTrend.select('residuals')
    .arrayProject([0])
    .arrayFlatten([['residuals']]);
  //print('harmonicTrendResiduals', harmonicTrendResiduals);
  
  // Pull out the bands to be visualized
  var sin = harmonicTrendCoefficients.select(ee.String(bandsForViz.get(0)));
  var cos = harmonicTrendCoefficients.select(ee.String(bandsForViz.get(1)));
  //print('sin', sin);
  //print('cos', cos);
  
  // Do some math to turn a first-order Fourier model into
  // hue, saturation, and value in the range[0,1].
  var amplitude = cos.hypot(sin).multiply(5).rename('amplitude');
  var phase = sin.atan2(cos).unitScale(-Math.PI, Math.PI).rename('phase');
  var val = harmonicImgColl.select(dependent).reduce('mean');
  //print('val', val);
  
  // Turn the HSV data into an RGB image and add it to the map.
  var seasonality = ee.Image.cat(phase, amplitude, val).hsvToRgb()
                            .rename(['seasonality_R',
                                     'seasonality_G',
                                     'seasonality_B']);
    
  // get SD of seasonality phase within sliding-window neighborhoods
  var phaseSD = phase.reduceNeighborhood({
    reducer: ee.Reducer.stdDev(),
    kernel: ee.Kernel.circle({radius: 100000, units: 'meters'})
  }).rename('phaseSD');
  //print('phaseSD', phaseSD);
  
                           
  //---------------------------------
  // Calculate the fitted vals and R2
  //---------------------------------
  
  var fittedHarmonic = harmonicImgColl.map(function(img) {
    return ee.Image(img).addBands(
      img.select(independents)
         .multiply(harmonicTrendCoefficients)
         .reduce('sum')
         .rename('fitted'));
  });
  //print('fittedHarmonic', fittedHarmonic);
  
  // Fix the time (in radians) of the maximum fitted value
  var tmax = fittedHarmonic.select(['fitted', 't'])
    .reduce(ee.Reducer.max(2))
    .rename('max_fitted', 't_max')
    .select('t_max')
    .mod(6.283185307179586);
    
  // Calculate the R2s
  var R2 = calcR2(fittedHarmonic, dependent);
  //print('R2', R2);
  
  // compose final output ImageCollection for the model
  var summary = seasonality.addBands(phase)
                           .addBands(amplitude)
                           .addBands(phaseSD)
                           .addBands(harmonicTrendCoefficients)
                           .addBands(harmonicTrendResiduals)
                           .addBands(tmax)
                           .addBands(R2)
                           // add max possible time_start val as property, to make it
                           // easy later to separate from original ImgColl it is merged with
                           .set({'system:time_start': 999999999999});
  //print('summary', summary);
  
  // merge the summary image and the rest of the ImageCollection
  var harmonicReg = ee.ImageCollection(summary)
    .merge(fittedHarmonic)
    // reproject to the input dataset's CRS
    //.map(function(img){return img.reproject(imgColl.first().projection())});
  return harmonicReg;
};

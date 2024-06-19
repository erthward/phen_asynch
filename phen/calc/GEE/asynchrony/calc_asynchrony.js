// TODO:

//  - devise my own home-grown function for calculating the regression stats (R2 mainly) as well

//  - inspect more deeply, to make sure that things make sense. Namely:
//    - Why do some scatterplots appear to have R2 values far below 1 at 
//      effectively 0 m distance (e.g. Point (-73.267, 9.805) at 2Km/px: R2 = 0.849 at 0.041 m)?


/*
  EXPORTS:
    calcAsynchrony(tsImgColl, maxNeighDist, minPctValidNeigh, proj)
      tsImgColl: The time-series ImageCollection to calculate asynchrony for
      maxNeighDist: The maximum distance (in meters) at which to include neighbors of a pixel;
                    used as the radius of a circular neighborhood kernel
      minPctValidNeigh: The minimum number of neighbors of a pixel that must be valid in order
                        for that pixel to be included in the output asynchrony Image.
      proj:             The CRS within which to run the calculation.
*/




//............................................................................................
// get R2s
//........

// function that takes a time-series ImageCollection, the max number of neighor cell-steps
// (the radius of the circular neighborhood kernel, expressed in pixels), and the CRS
// to work in, and returns an ImageCollection of R2s, where each Image contains the R2s
// of the correlations between the pixels and their neighbors at a certain position
// (such that the size of the ImageCollection equals
//  the number of neighbors being used to calculate asynchrony)
var calcNeighborTimeSeriesR2s = function(tsImgColl, maxNeighSteps, proj){

  // convert time series to an array per pixel
  // (NOTE: I checked and it becomes a 101 x 1 array image, i.e. n_ts_images x n_ts_bands)
  var tsArray = tsImgColl.toArray();
  
  // collect neighborhood into bands
  // (NOTE: for SIF with maxNeighSteps=?, produces an Image of 317 bands,
  //  each one 101 x 1 array image,
  //  i.e. maxNeighSteps' resulting-neighbor-count-number of bands,
  //  with each band being an n_ts_images x n_ts_bands array)
  var tsNeighBands = tsArray
    // then grab all the neighbors of each pixels into band-values at that pixel
    .neighborhoodToBands(ee.Kernel.circle(maxNeighSteps, 'pixels'));
  
  // get the time series for the focal pixels
  var tsFocal = tsNeighBands.select('array_0_0');
  
  // convert the multiband Image into a multi-image ImageCollection
  var tsBandNames = tsNeighBands.bandNames();
  var tsNeighImgColl = ee.ImageCollection.fromImages(tsBandNames.map(function(name){
    name = ee.String(name);
    var img = tsNeighBands.select(name).toFloat()
      .set('origBandName', name);
    return img;
  }));
  
  // function that takes a neighbor-array image, concatenates its array with the focal-cell
  // array-image, then reduces the concatenated array image to get the R2 between the neigh's
  // and focal cell's time series
  var calcCorrelation = function(neighImg){
    // concatenate the two arrays into a n_ts_images x 2 array-image
    var focalAndThisNeighCatArray = neighImg.arrayCat(tsFocal, 1);
    // array-reduce that array image using pearsonsCorrelation,
    // (reduce across the 0 axis, using the 1 axis as the input fields),
    // then get the Pearson's corr coeff (item [0,0] in the output array) and square it to get R2
    var focalAndThisNeighR2Img = focalAndThisNeighCatArray
      .arrayReduce(ee.Reducer.pearsonsCorrelation(), [0], 1).arrayGet([0,0]).pow(2)
      .set('origBandName', ee.String(ee.Image(neighImg.get('origBandName'))));
    return focalAndThisNeighR2Img;
  };
  
  // map the neighbor-correlation function across the IC, allowing nulls to be dropped
  var R2sNeedReproj = tsNeighImgColl.map(calcCorrelation, true);
  
  // reproject to make sure that calculations happen in the input data's CRS
  var R2s = R2sNeedReproj.map(function(img){return img.reproject(proj)});
  
  return R2s;
};




//...........................................................................................
// get dists
//..........

// function that takes the max number of neighor cell-steps
// (the radius of the circular neighborhood kernel, expressed in pixels), and the CRS
// to work in, and returns an ImageCollection of distances (in meters),
// where each Image contains the dists between the pixels and their neighbors at a certain position
// (such that the size of the ImageCollection equals
//  the number of neighbors being used to calculate asynchrony)
var calcNeighborDists = function(maxNeighSteps, proj){

  // grab pixelCoordinates as lon, lat
  var coordArray = ee.Image.pixelLonLat()
    // convert to radians
    .divide(180).multiply(3.141592654)
    // reduce to an array per pixel
    .toArray();
    
  // collect neighbors into bands
  var coordNeighBands = coordArray.neighborhoodToBands(ee.Kernel.circle(maxNeighSteps, 'pixels'));
  
  // get the coordinates for the focal pixels, converted to radians
  var coordsFocal = coordNeighBands.select('array_0_0');
  
  // convert the multiband Image into a multi-image ImageCollection
  var coordBandNames = coordNeighBands.bandNames();
  var coordNeighImgColl = ee.ImageCollection.fromImages(coordBandNames.map(function(name){
    name = ee.String(name);
    var img = coordNeighBands.select(name).toFloat()
      .set('origBandName', name);
    return img;
  }));
  
  // grab the lats as their own separate ImageCollections
  // (to be combined with latDiffs and lonDiffs later)
  var latNeighImgColl = coordNeighImgColl.map(function(img){
    return img.arrayGet(1)});
  
  // map the distance-from-focal-cell calculation across all images
  // (i.e. neighbors) in the collection, using haversine formula to get answer in km
  var diffs = coordNeighImgColl.map(function(img){
    return ee.Image(img).subtract(coordsFocal)});
  var lonDiffs = diffs.map(function(img){
    return img.arrayGet(0)});
  var latDiffs = diffs.map(function(img){
    return img.arrayGet(1)});
  
  // combine all of the coordinate-derived variables needed for the Haversine
  // formula into a single ImageCollection
  var combinedCoordVars = latDiffs.combine(latNeighImgColl).combine(lonDiffs)
    .map(function(img){
      var origBandName = ee.String(img.bandNames().get(0));
      var combinedImg = img.rename('latDiff', 'lat', 'lonDiff')
        .set('origBandName', origBandName);
      return combinedImg});
  
  // calculate distances in meters, using haversine formula
  var distsNeedReproj = combinedCoordVars.map(function(img){
    return ee.Image(ee.Image(img).select('latDiff').divide(2).sin().pow(2)
        .add(ee.Image(coordsFocal.arrayGet(0)).cos()
          .multiply(img.select('lat').cos())
          .multiply(img.select('lonDiff').divide(2).sin().pow(2))))
      // NOTE: 6371008 m is a good approximation of any of the three common methods for
      //       representing the Earth by its 'mean' radius
      // NOTE: would technically need to ensure that the component inside the sqrt does not
      //       exceed 1 (to avoid getting non-real answers), but I don't actually need to
      //       worry about that here given that that only actually happens for points
      //       approaching antipodes, which is far out  side what I'll be using here
      .sqrt().asin().multiply(2).multiply(6371008)
      .rename(ee.String(ee.Image(img).get('origBandName')))
      .set('origBandName', ee.String(ee.Image(img).get('origBandName')))});
  
  // reproject to make sure that calculations happen in the input data's CRS
  var dists = distsNeedReproj.map(function(img){return img.reproject(proj)});
 
  return dists;
};



//............................................................................................
// correlate the two resulting ImageCollections
//.............................................


// NOTE: The linearFit reducer won't allow me to omit (i.e. fix) the intercept (i.e. 'offset'),
// and the linearRegression reducer outputs too much stuff, so won't work with an array image.
// So I define and run my own linear-regression function instead
// (one that allows me to feed in the value at which to fix the intercept)
var calcOLSRegressionWithFixedIntercept = function(x, y, intercept){
  var X = ee.Image(x);
  var Y = ee.Image(y);
  var INTERCEPT = ee.Number(intercept);
  var betas = X
    .matrixTranspose(0,1)
    .matrixMultiply(X)
    .matrixInverse()
    .matrixMultiply(X.matrixTranspose())
    .matrixMultiply(Y.subtract(INTERCEPT))
    .arrayGet([0,0]);
  return betas;
};


// function that takes a time-series ImageCollection, the maximum neighbor-distance to include
// in the each pixel's neighborhood (used as the radius of a circular kernel),
// the minimum number of valid neighbors that a pixel needs in order to be included in the output,
// and the CRS to operate it, and returns an Image of calculate asynchrony values,
// with pixels less than minPctValidNeigh valid neighbors being masked out
exports.calcAsynchrony = function(tsImgColl, maxNeighDist, minPctValidNeigh, proj){
  
  // get the max neighbor-distance in pixels ('nNeigh')
  // using the max neighbor-distance in meters ('maxNeighDist')
  var maxNeighSteps = ee.Number(maxNeighDist).divide(proj.nominalScale()).ceil();
  
  // calculate the neigh-R2s ImageCollection  
  var R2s = calcNeighborTimeSeriesR2s(tsImgColl, maxNeighSteps, proj);
  
  // calculate the min number of valid neighbor pixels ('minValidNeigh'),
  // using the min percent of valid neighbors ('minPctValidNeigh')
  var minValidNeigh = ee.Number(minPctValidNeigh).multiply(ee.Number(R2s.size())).ceil();
  
  // calculate the neigh-dists ImageCollection
  var dists = calcNeighborDists(maxNeighSteps, proj);
  
  // sort the R2s and dists so that they're in identical order
  var R2sSorted = R2s.sort('origBandName')
    // set the bandname to an identical, non-descript bandname
    // (to avoid downstream problem arising from non-homogeneously named ImageCollection)
    .map(function(img){
      return img.rename('identical_bandname')});
  var distsSorted = dists.sort('origBandName')
    // set the bandname to an identical, non-descript bandname
    // (to avoid downstream problem arising from non-homogeneously named ImageCollection)
    .map(function(img){
      return img.rename('identical_bandname')});
      
  // use the R2s ImageCollection to grab a mask based on the minimum number
  // of valid neighbors a pixel must have ('minValidNeigh'),
  // by counting non-null R2s at each pixel
  // (NOTE: the mask will be applied after the asynchrony calculation)
  var minValidNeighMask = R2sSorted.count().gte(ee.Number(minValidNeigh));
  
  // combine the sorted R2s and dists ImageCollections
  // into a single, 2-band ImageCollection,
  var distsAndR2sSorted = distsSorted.combine(R2sSorted).map(function(img){
    return img.rename('dist', 'R2')});

  // then use the combined ImageCollection to mask the dists wherever they
  // will have no paired R2 value because of data dropout
  // NOTE: I THINK THIS ONLY HAPPENS BECAUSE I PUT THE FITTED VALUES INTO THE ORIGINAL
  // SIF DATASET?
  var distsSortedMasked = distsAndR2sSorted.map(function(img){
    var R2 = img.select('R2');
    return img.select('dist').mask(R2.mask())});
  
  // call the fixed-intercept regression function on the dists and R2s
  var asynchrony = calcOLSRegressionWithFixedIntercept(distsSortedMasked.toArray(), 
                                                        R2sSorted.toArray(), 
                                                        1)
    // mask the pixels that didn't have >= min number of valid neighbors
    .updateMask(minValidNeighMask)
    // take the absolute value of the slope as our metric
    .abs()
    // and rename the band
    .rename('asynch');
  
  return asynchrony;
};



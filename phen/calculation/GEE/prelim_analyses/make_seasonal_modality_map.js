//--------
// imports
//--------

var hReg = require('users/drewhart/seasonality/:fft/harmonic_regression.js');
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
var cbar = require('users/drewhart/seasonality/:viz/make_colorbar.js');

//---------
// palettes
//---------

// load palettes
var palettes = require('users/gena/packages:palettes'); 


//----------
// NIRv data
//----------

var NIRvP = gDat.getNIRvPData();
//print('NIRv', NIRv.first());


//---------
// SIF data
//---------

// read in SIF data and reformat its time data
var SIF = gDat.getSIFData();
//print('SIF', SIF);



//----------------------
// run regressions: NIRv
//----------------------

// calculate regressions for annual seasonality
var dependent = 'NIRvP';
var harmonics = ee.List([1]);
var bandsForViz = ee.List(['sin_1', 'cos_1']);
var annNIRvP = ee.ImageCollection(hReg.calcHarmonicRegression(NIRvP.limit(1000), dependent,
                                                         harmonics, bandsForViz));

// calculate regressions for semiannual seasonality
var dependent = 'NIRvP';
var harmonics = ee.List([2]);
var bandsForViz = ee.List(['sin_2', 'cos_2']);
var semiNIRvP = ee.ImageCollection(hReg.calcHarmonicRegression(NIRvP.limit(1000), dependent,
                                                          harmonics, bandsForViz));

// compare both the regressions' R2 images
var resultNIRvP = ee.Image(semiNIRvP.first().select('R2'))
  .gte(ee.Image(annNIRvP.first().select('R2')))
  // add 1, so that the resulting image is 1 for annual (1 cycle per year),
  // 2 for semiannual (2 cycles per year);
  .add(ee.Number(1));
var diffNIRvP = ee.Image(semiNIRvP.first().select('R2'))
  .subtract(ee.Image(annNIRvP.first().select('R2')))
  .abs();


/*
//-------
// export
//-------

//var geom = ee.Geometry.Rectangle([[-180, -90], [180, 90]]);
Export.image.toAsset({
  image: result,
  description: 'seasonalModeNIRv',
  assetId: 'seasonalModeNIRv',
  scale: NIRv.first().projection().nominalScale().getInfo(),
//  region: geom,
  pyramidingPolicy: {
    'R2': 'mode',
  },
  maxPixels: 4000000000
});
*/


//---------------------
// run regressions: SIF
//---------------------

// calculate regressions for annual seasonality
var dependent = 'b1';
var harmonics = ee.List([1]);
var bandsForViz = ee.List(['sin_1', 'cos_1']);
var ann = ee.ImageCollection(hReg.calcHarmonicRegression(SIF, dependent,
                                                         harmonics, bandsForViz));

// calculate regressions for semiannual seasonality
var dependent = 'b1';
var harmonics = ee.List([2]);
var bandsForViz = ee.List(['sin_2', 'cos_2']);
var semi = ee.ImageCollection(hReg.calcHarmonicRegression(SIF, dependent,
                                                          harmonics, bandsForViz));

// compare both the regressions' R2 images
var result = ee.Image(semi.first().select('R2'))
  .gte(ee.Image(ann.first().select('R2')))
  // add 1, so that the resulting image is 1 for annual (1 cycle per year),
  // 2 for semiannual (2 cycles per year);
  .add(ee.Number(1));

var diff = ee.Image(semi.first().select('R2'))
  .subtract(ee.Image(ann.first().select('R2')))
  .abs();

/*
//-------
// export
//-------

//var geom = ee.Geometry.Rectangle([[-180, -90], [180, 90]]);
Export.image.toAsset({
  image: result,
  description: 'seasonalModeSIF',
  assetId: 'seasonalModeSIF',
  scale: SIF.first().projection().nominalScale().getInfo(),
//  region: geom,
  pyramidingPolicy: {
    'R2': 'mode',
  },
  maxPixels: 4000000000
});

*/



//-------
// map it
//-------
Map.addLayer(resultNIRvP, {palette: ['green', 'orange'], min:1, max:2, opacity: 0.8}, 'result (calculated live)');
Map.addLayer(diffNIRvP, {palette: palettes.cmocean.Turbid[7].reverse(), min:0, max:1, opacity: 0.8}, 'diff (calculated live)');
Map.addLayer(SIF, {opacity: 0.05}, 'SIF');
Map.addLayer(NIRvP.limit(1000), {opacity: 0.05}, 'NIRvP');
//var SIFresult = ee.Image('users/drewhart/seasonalModeSIF');
//print(SIFresult);
//var NIRvresult = ee.Image('users/drewhart/seasonalModeNIRv');
//print(NIRvresult);

// default to terrain
Map.setOptions('TERRAIN');
//Map.addLayer(SIFresult, {palette: palettes.kovesi.isoluminant_cgo_70_c39[7],
//                      min: 1, max: 2, opacity: 0.8},
//  'SIF seasonal modality (blue: annual, pink: semiannual)');
//Map.addLayer(NIRvresult, {palette: palettes.matplotlib.viridis[7],
//                      min: 1, max: 2, opacity: 0.8},
//  'NIRv seasonal modality (navy: annual, yellow: semiannual)');



// make map of time of year of max fitted value
// (sort of a gross approx. of 'phase' that applies to both annual
// and semiannual sites
var dependent = 'NIRvP';
var harmonics = ee.List([1, 2]);
var bandsForViz = ee.List(['sin_1', 'cos_1']);
var full = ee.ImageCollection(hReg.calcHarmonicRegression(NIRvP, dependent,
                                                         harmonics, bandsForViz));
print('result', full.limit(10));
Map.addLayer(full.first().select('t_max'), 
  {palette: palettes.kovesi.cyclic_mygbm_30_95_c78_s25[7], opacity: 0.8,
   min: 0, max: 6.283185307179586},
  't max');
Map.add(cbar.makeLegend('annual-frequency phase', 'Jan 1', 'Jun 1', 'Dec 31', 
  palettes.kovesi.cyclic_mygbm_30_95_c78_s25[7], 'top-right'));
 
/* 
//-------
// export
//-------

//var geom = ee.Geometry.Rectangle([[-180, -90], [180, 90]]);
Export.image.toAsset({
  image: full.first().select('t_max'),
  description: 'fittedpeakNIRv',
  assetId: 'fittedpeakNIRv',
  scale: NIRv.first().projection().nominalScale().getInfo(),
//  region: geom,
  pyramidingPolicy: {
    't_max': 'mean',
  },
  maxPixels: 4000000000
});
*/
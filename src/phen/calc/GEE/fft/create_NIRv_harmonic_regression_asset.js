// script to produce the NIRv harmonic regression and save as an asset 


//------------
// LOAD PARAMS
//------------
var params = require('users/drewhart/seasonality/:params.js');


//------------------------
// LOAD REQUIRED FUNCTIONS
//------------------------
// data-loading
var gDat = require('users/drewhart/seasonality/:io/get_data.js');

// harmonic regression
var hReg = require('users/drewhart/seasonality/:fft/run_harmonic_regression.js');


//----------
// LOAD DATA
//----------
var NIRv = gDat.getNIRvData(params.maskLowBRDFAlbQual,
                            params.maxBRDFAlbQualVal,
                            params.maskWaterSeparately,
                            params.NIRvNYears,
                            params.NIRvEndYear,
                            params.NIRvDayStep,
                            params.NIRvClampToMinPos,
                            params.NIRvUnmaskToMinPos,
                            params.NIRvMaskLteZero
                           );
if (params.verbose){
  print('NIRvDayStep', params.NIRvDayStep);
  print('NIRv', NIRv);
}
var dependent = 'NIRv';
var proj = NIRv.first().projection();
var scale = NIRv.first().projection().nominalScale();


//---------------------
// CALCULATE REGRESSION
//---------------------
// get an ImageCollection where the first Image is the regression summary (coeffs and R2)
// and the following 365 images are a one-year daily time series of fitted values derived from that regression
var reg = ee.ImageCollection(hReg.calcHarmonicRegression(NIRv, 
                                                         dependent,
                                                         params.harmonics, 
                                                         params.bandsForViz));


//----------
// EXPORT IT
//----------
var scale = NIRv
  .first()
  .projection()
  .nominalScale()
  .getInfo();
var assetName = 'NIRv_harmonic_regression';
Export.image.toAsset({
  // NOTE: export just the first image with coefficients and R2
  image: reg.first(),
  description: assetName,
  assetId: assetName,
  scale: NIRv.first().projection().nominalScale().getInfo(),
  region: params.roi,
});

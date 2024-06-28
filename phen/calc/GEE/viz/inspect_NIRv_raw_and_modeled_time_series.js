var params = require('users/drewhart/seasonality/:params.js');
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
var hReg = require('users/drewhart/seasonality/:fft/run_harmonic_regression.js');

// use just 4 years of data to plot model results
var NIRv = gDat.getNIRvData(params.maskLowBRDFAlbQual,
                            3,
                            false,
                            4,
                            2021,
                            7,
                            params.NIRvClampToMinPos,
                            params.NIRvUnmaskToMinPos,
                            params.NIRvMaskLteZero
                           );
var dependent = 'NIRv';
print('NIRv', NIRv);
Map.addLayer(NIRv.select('NIRv'), {}, 'NIRv', false);

var LSP = ee.ImageCollection(hReg.calcHarmonicRegression(NIRv, 
                                                         dependent,
                                                         params.harmonics, 
                                                         params.bandsForViz));                                     
print('LSP', LSP);
// Map.addLayer(LSP.select('fitted'), {}, 'LSP reg', false);
Map.addLayer(ee.ImageCollection(LSP.toList(365*4,1)).select(['NIRv', 'fitted']), {}, 'LSP fit', false);
Map.addLayer(LSP.first(), {}, 'LSP reg results', false);

Map.addLayer(LSP.first().select('R2'), {}, 'R2');
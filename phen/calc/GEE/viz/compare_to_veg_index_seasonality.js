var phen =ee.ImageCollection("MODIS/006/MCD12Q2");
print('phen', phen);

var palettes = require('users/gena/packages:palettes');
var cbar = require('users/drewhart/seasonality/:viz/make_colorbar.js');

Map.addLayer(phen
  .select('NumCycles')
  .reduce(ee.Reducer.mean()),
  {palette: ['green', 'orange'], min:1, max:2}, 'n cycles');
  
Map.addLayer(phen
  .select('Peak_1')
  .map(function(img){return img.mod(365)})
  .reduce(ee.Reducer.mean()),
  {palette: palettes.kovesi.cyclic_mygbm_30_95_c78[7], min:0, max:365}, 'peak 1');

Map.add(cbar.makeLegend('phase', 'Jan 1', 'Jul 1', 'Dec 31',
          palettes.kovesi.cyclic_mygbm_30_95_c78_s25[7], 'bottom-center'));
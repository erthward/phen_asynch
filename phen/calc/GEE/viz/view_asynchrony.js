///////////////////////////////////////////////
// plot asynchrony
var asynch = ee.Image("users/drewhart/global_phenol_asynch_SIF");
print('global asynchrony, ANN SIF', asynch);

var asynch_NIRvP = ee.Image("users/drewhart/NA_SA_agg_asynch");
print('asynchrony, NA and SA, NIRvP', asynch_NIRvP);


//TODO
// 1. figure out why the ingested asset isn't properly masked
// 2. figure out how to set band names correctly




// take care of masking issue
var asynch = asynch.updateMask(asynch.abs().gte(0));

var cbar = require('users/drewhart/seasonality/:viz/make_colorbar.js');
var palettes = require('users/gena/packages:palettes'); 

var maxVal = 0.000018;

var custom_palette = false;

Map.setOptions('TERRAIN');
// get the palettes
if (custom_palette){
  var asynch_pal = ['#000000', '#171616', '#362f2f', '#523a3a', 
               '#7d3e3e', '#9e3737', '#c42727', '#f52525'];
} else {
  var asynch_pal = palettes.crameri.imola[50];
}
var r2_pal = palettes.crameri.hawaii[50];
var n_pal = palettes.crameri.acton[50];

// add the black background
Map.addLayer(ee.Image.constant(0), {min:0, max:1}, 'black background', true);

// add the asynch layers (R2-based and Euclidean distance-based)
Map.addLayer(asynch.select('b1'), {min:0, max: 0.0000076, opacity:0.8, palette: asynch_pal},
             'asynchrony (R2-based)', false);

Map.addLayer(asynch.select('b3'), {min:0.000027, max: 0.00027, opacity:0.8, palette: asynch_pal }, 
             'asynchrony (Euclidean distance-based), SIF', true);

Map.addLayer(asynch_NIRvP.select('b3'), {min:-0.000046, max: 0.00013, opacity:0.8, palette: asynch_pal }, 
             'asynchrony (Euclidean distance-based), NIRvP', false);             
             
// add the R2s for the 2 asynch layers
Map.addLayer(asynch.select('b2'), {min:0, max: 1, opacity:0.8,
                                   palette: r2_pal},
             'R2 (R2-based asynchrony)', false);
Map.addLayer(asynch.select('b4'), {min:0, max: 1, opacity:0.8,
                                   palette: r2_pal},
             'R2 (Euclidean distance-based asynchrony)', false);
             
// add the asynch sample sizes
Map.addLayer(asynch_NIRvP.select('b5'), {min:1200, max: 4000, opacity:0.8,
                                   palette: n_pal},
             'sample size (n pixels use to calculate asynch', false);
             
// add the colorbar
Map.add(cbar.makeLegend('asynchrony', 'low', 'med', 'high', asynch_pal, 'bottom-center'));
//Map.add(cbar.makeLegend('R2', '<=0', '0.5', '1', r2_pal, 'bottom-left'));
//Map.add(cbar.makeLegend('sample size', '1200', 'med', '4000', n_pal, 'bottom-right'));

///////////////////////////////////////////////
// add landcover
var modis = ee.ImageCollection("MODIS/006/MCD12Q1")
  .select('LC_Type1')
  .filterDate('2014-01-01', '2019-01-01');
var LC14 = modis.first();
var LC14Viz = {
  min: 1.0,
  max: 17.0,
  opacity: 0.5,
  palette: [
    '05450a', '086a10', '54a708', '78d203', '009900', 'c6b044', 'dcd159',
    'dade48', 'fbff13', 'b6ff05', '27ff87', 'c24f44', 'a5a5a5', 'ff6d4c',
    '69fff8', 'f9ffa4', '1c0dff'
  ],
};
Map.addLayer(LC14, LC14Viz, 'LC', false);


var modis_lc = ee.ImageCollection("MODIS/006/MCD12Q1")
  .filterDate('2018-12-31', '2019-01-02');
print('MODIS LC', modis_lc);
Map.addLayer(modis_lc.select('LC_Type1'), {min:0, max:11, palette: palettes.colorbrewer.Set3[12]},
             'land cover type', false);
Map.addLayer(ee.Image(modis_lc.select('LC_Type1').first())
  .neq(ee.Image.constant(3))
  .multiply(ee.Image(modis_lc.select('LC_Type1').first()).neq(ee.Image.constant(4))),
  {min:0, max:1, palette: ['black', 'green']}, 'not deciduous forest', false);


// add biodiversity hotspots
var spots = ee.FeatureCollection("users/drewhart/hotspots");
print('biodiversity hotspots', spots);

var empty = ee.Image().byte();
var outline = empty.paint({
  featureCollection: spots,
  color: 1,
  width: 1
});
Map.addLayer(outline, {palette: '#f7f7f7'}, 'biodiversity hotspots', false);

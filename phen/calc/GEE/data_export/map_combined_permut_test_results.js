// read in the combined results
var res = ee.Image('users/drewhart/testRes_0-299_combo');

var ntests = ee.Number(300);

// alpha, to determine significance (for mapping only)
var alpha = ee.Number(0.01);
var threshold = alpha.multiply(ntests).ceil().getInfo();

// map it
var palettes = require('users/gena/packages:palettes'); 
Map.addLayer(res, {min: 0, max: threshold,
                   palette: palettes.crameri.roma[10].reverse(), opacity: 0.8},
             'combined results');

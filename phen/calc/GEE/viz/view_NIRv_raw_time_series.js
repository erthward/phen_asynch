var gDat = require('users/drewhart/seasonality/:io/get_data.js');
var NIRv = gDat.getNIRvData(50, 10);
Map.addLayer(NIRv.map(function(img){return img.select('NIRv')}));
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
var sif = gDat.getSIFData(100);
var sif_missing = sif.reduce(ee.Reducer.count()).unmask().lte(ee.Number(0.1).multiply(sif.size()));

var nirv = gDat.getNIRvData(100);
var nirv_not_missing = nirv.reduce(ee.Reducer.count()).unmask().gte(ee.Number(0.5).multiply(nirv.size()));
var sif_missing_nirv_not = sif_missing.and(nirv_not_missing
  .setDefaultProjection({crs: nirv.first().projection(), scale: nirv.first().projection().nominalScale()}))
  .reduceResolution(ee.Reducer.mode())
  .reproject({crs:sif.first().projection(), scale:sif.first().projection().nominalScale()});
Map.addLayer(sif_missing, {}, 'sif missing');
Map.addLayer(nirv_not_missing, {}, 'nirv not');N
Map.addLayer(sif_missing_nirv_not, {}, 'both');

Map.setZoom(0);
var roi = ee.Geometry.Polygon(
        [[[-163.34629019414842, 75.64317543582969],
          [-163.34629019414842, -63.28634219669422],
          [186.80995980585158, -63.28634219669422],
          [186.80995980585158, 75.64317543582969]]], null, false);
Export.image.toAsset({image: sif_missing_nirv_not,
                      description: 'SIF_missing_NIRv_not',
                      assetId: 'users/drewhart/SIF_missing_NIRv_not',
                      pyramidingPolicy: 'mean',
                      maxPixels: 4000000000,
                      region: roi,
                      scale: sif.first().projection().nominalScale().getInfo(),
                     });

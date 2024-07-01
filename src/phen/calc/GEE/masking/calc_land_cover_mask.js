// get necessary parameters
var params = require('users/drewhart/seasonality/:params.js');
var end_yr = ee.Number(params.NIRvEndYear);
var start_yr = end_yr.subtract(ee.Number(params.NIRvNYears));
var date_ending = ee.String('-01-01T00:00');
var start_datestr = ee.String(start_yr).cat(date_ending);
var end_datestr = ee.String(end_yr).cat(date_ending);
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
var SIF = gDat.getSIFData(false);
var SIFProj = SIF.first().projection();
var SIFScale = SIFProj.nominalScale();

// create a remapping of all LC classes,
// to define those to be included in both LSP and asynch analyses (new classes 1 and 2)
// and those to be included only in LSP analyses (new class 1)
var validLCClasses = ee.List([1,2,3,4,5,6,7,8,9,10,11,12,14]);
var remapClasses = ee.List([1,1,1,1,1,1,1,1,1,1,1,2,2]);

// load land cover dataset
// NOTE: replaced previously used MODIS-res product (MCD12Q1)
// with pre-aggregated majority LC product recommended by Reviewer 1 (MCD12C1)
var lc = ee.ImageCollection("MODIS/061/MCD12C1")
  .filterDate(start_datestr, end_datestr)
  .select('Majority_Land_Cover_Type_1');

// create a separate, simple water mask asset
// (to be used below and to be used for plotting supp figs of all mask maps)
var waterMask = lc
  .map(function(img){return img
    .reduceResolution({
        reducer: ee.Reducer.mode(),
        maxPixels: 2000})
      .reproject({crs: SIFProj, scale: SIFScale})
    .gt(0)
    .rename('waterMask')})
  .reduce(ee.Reducer.allNonZero());

// remap pixels to the classes defined above
var lcReclass = lc
  .map(function(img){return img
    .reduceResolution({
        reducer: ee.Reducer.mode(),
        maxPixels: 2000})
      .reproject({crs: SIFProj, scale: SIFScale})
    // mask out all water off the bat
    .updateMask(waterMask)
    .remap({from: validLCClasses,
            to: remapClasses,
            defaultValue: 0
           })
    .rename('lcMask')});

// reduce time series, flagging pixels that are always 'natural'
// (even if they change class, since a lot of hard-to-class regions do so meaninglessly,
// e.g., western Chaco, woodland and xeric regions of southern Africa)
// and/or always agricultural (the latter of which will be dropped only from asynch analyses,
// since they could easily show anthropogenic asynchrony in veg dynamics)
var alwaysNat = lcReclass.map(function(img){return img.eq(1)}).reduce(ee.Reducer.allNonZero());
var alwaysAg = lcReclass.map(function(img){return img.eq(2)}).reduce(ee.Reducer.allNonZero());

// combine into a single mask map, in which 1s are either always 'natural' or always agricultural
// and 2s are strictly always natural (such that all non-zeros will be used in LSP analyses
// but only 2s will be used in asynch analyses)
var lcMask = alwaysNat.or(alwaysAg).add(alwaysNat);

// map the results, if requested
if (params.map_intermediates){
  Map.addLayer(lc, {}, 'original MCD12C1 land cover', false);
  Map.addLayer(lcMask.gt(0).selfMask(), {palette:['black', 'yellow'], opacity:0.7}, 'lcMask: default', true);
  Map.addLayer(lcMask.gt(1).selfMask(), {palette:['black', 'red'], opacity:0.7}, 'lcMask: strict', true);

  }

// export masking layers
Export.image.toAsset({image: waterMask,
                      description: 'LSP_waterMask',
                      assetId: 'users/drewhart/LSP_waterMask',
                      region: params.roi,
                      pyramidingPolicy: 'mode',
                      maxPixels: params.maxPixels,
                      scale: lc.first().projection().nominalScale().getInfo()});
Export.image.toAsset({image: lcMask,
                      description: 'LSP_lcMask',
                      assetId: 'users/drewhart/LSP_lcMask',
                      region: params.roi,
                      pyramidingPolicy: 'mode',
                      maxPixels: params.maxPixels,
                      scale: lc.first().projection().nominalScale().getInfo()});

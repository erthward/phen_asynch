var dataset= 'NIRv';
//var dataset = 'SIF';

var ts = ee.Image('users/drewhart/NIRv_ann_fit_vals');

var bandNames = ts.bandNames();
var bandNums = ee.List.sequence(1, bandNames.size());
var tsImgColl = ee.ImageCollection(bandNames.map(function(name){
  return ts.select([name]).rename('fitted');
}));

// GLOBAL MP4 VIDEO EXPORT CODE

/// CREATE RGB VISUALIZATION IMAGES FOR VIDEO EXPORT ///
// creates a continumn of colors for seasonality animation
var rgbVis = tsImgColl.map(function(img) {
  return img.visualize(
        {
        min: 0,
        max: 1,
        palette: [
          'B6DC1A','A8CA19','99B718', '74A901','66A000', 
          '529400','3E8601', '207401', '056201', '004C00', 
          '023B01'
        ]}
        )});


// Need to have 3-band imagery for the video.
var rgbVis = rgbVis.select(['vis-red', 'vis-green', 'vis-blue']);
print(rgbVis);

// export params
var globalROI = 
    ee.Geometry.Polygon(
        [[[-180, -90],
          [180, -90],
          [180, 90],
          [-180, 90],
          [-180, -90]]], null, false);


/// Export video to Drive ///
/* NOTE: CANNOT HAVE PARAMETER VALUES IN BOTH 'scale' & 'dimensions'
    CAN ONLY HAVE ONE OR THE OTHER */
Export.video.toDrive({
  folder: 'phen_outputs_from_GEE', 
  // use 'collection' image collection for full resolution
  collection: rgbVis, 
  description: ee.String(dataset).cat(ee.String('_ann_fit_vals_vid')).getInfo(),
  // 'scale' = the resolution of the respective image collection
  //scale: collection.first().projection().nominalScale().getInfo(), // meters per pixel 
  framesPerSecond: 7, 
  maxPixels: 8e9,
  maxFrames: 1000,
  dimensions: 800,
  crs: rgbVis.first().projection().crs().getInfo(),
  region: globalROI 
}); 



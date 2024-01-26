/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var image = ee.Image("USGS/GMTED2010"),
    image2 = ee.Image("users/nguyenbui/vrm_1KMmn_GMTEDmd"),
    image4 = ee.Image("USGS/SRTMGL1_003"),
    image5 = ee.Image("users/nguyenbui/VRM10km");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
//Map.addLayer(image5)
//Testing out Terrain Ruggedness Index and Topographic Positionality Index
var dataset = ee.Image('USGS/GMTED2010');
var elevation = dataset.select('be75');
//Masking water
//var minOccurr = 80;
//var waterMask1 = ee.Image("JRC/GSW1_0/GlobalSurfaceWater")
//    .select('occurrence')
//    .gt(minOccurr)
//    .unmask(0)
//    .not();
//var waterMask2 = ee.ImageCollection("MODIS/006/MOD44W")
//  .select('water_mask')
//  .reduce(ee.Reducer.sum())
//  .gt(0)
//  .not();
// combine the two water masks
//var waterMask = waterMask1.and(waterMask2);
// mask out water in all images in the collection
var data = elevation;
var elevationVis = {
  min: -100.0,
  max: 6500.0,
  gamma: 3.5,
};
Map.setCenter(-122, 20, 8);
//Map.addLayer(data, elevationVis, 'Elevation');

//var vrm1km = image2.updateMask(waterMask)
//Map.addLayer(vrm1km, {min:0, max:0.1}, '1km VRM')

var slope = ee.Terrain.slope(data)
//Map.addLayer(slope, {min:0, max:90}, 'Slope');

var aspect = ee.Terrain.aspect(data)
//Map.addLayer(aspect, {min:0, max:360}, 'Aspect')

var z = slope.divide(180).multiply(Math.PI).cos()
var xy = slope.divide(180).multiply(Math.PI).sin()
var x = aspect.divide(180).multiply(Math.PI).sin().multiply(xy)
var y = aspect.divide(180).multiply(Math.PI).cos().multiply(xy)

var xyz = ee.Image.cat([x,y,z]).select(
    ['aspect', 'aspect_1', 'slope'], // old names
    ['x', 'y', 'z']               // new names
);
Map.addLayer(xyz, {min:-1, max:1}, 'xyz')
var projxyz = xyz.projection();
//var scalexyz = xyz.projection.scale();
//var txyz = xyz.projection.transform
//print(scalexyz)

//Input kernel size
var size = 50000

var sumxyz = xyz.reduceNeighborhood({
  reducer: ee.Reducer.sum().combine({
  reducer2: ee.Reducer.count(),
  sharedInputs: true}),
  kernel: ee.Kernel.circle(size,'meters',false)
}).reproject(projxyz)

//Map.addLayer(sumxyz, {min:-1, max:1}, 'sumxyz')

var ssxyz = sumxyz.select('x_sum').pow(2).add(sumxyz.select('y_sum').pow(2)).add(sumxyz.select('z_sum').pow(2))
//Map.addLayer(ssxyz, {min:-1, max:1}, 'ssxyz')

var r = ssxyz.sqrt()
//Map.addLayer(r, {min:-1, max:1}, 'r')

var vrm = r.divide(sumxyz.select('x_count')).multiply(-1).add(1)
Map.addLayer(vrm, {min:0, max:0.1}, 'vrm')
var scale = vrm.projection().nominalScale();
var roi =
    ee.Geometry.Polygon(
        [[[-165, 60],
          [-165, -60],
          [180, -60],
          [180, 60],
          [-165,60]]], null, false);

Export.image.toAsset({
  image:vrm,
  description:"VRM50km",
  scale:scale.getInfo(),
  maxPixels:100000000000,
  region:roi,
})

//some values are negative and over 1. Identify a way to identify all the places that are larger than 1 and less than 0, so i can fix the code


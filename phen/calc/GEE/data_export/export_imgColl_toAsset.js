
/*
FUNCTION: exportImageCollectionToAsset
- Adapted from:
    * Author: Rodrigo E. Principe
    * License: Apache 2.0
    * email: fitoprincipe82@gmail.com
    https://gis.stackexchange.com/questions/407146/export-imagecollection-to-asset
- DESCRIPTION: Exports a series of images for every n days from an 
    ImageCollection as a GEE Asset.
- PARAMETERS:
    * collection: (ImageCollection) an image collection
    * assetId: (string) the name you choose for the exported Assets
    * rio: (Geometry.LinearRing|Geometry.Polygon|String) region of interest 
    * scale: (number) the scale or meters per pixel of the ImageCollection 
    * crs: (string) the CRS of the ImageCollection 
    * maxPixels: (number) Restrict the number of pixels in the export. By 
        default, you will see an error if the export exceeds 1e8 pixels
    * start: (number) starting index
    * stepSize: (number) step size
- PARAMETER DEFAULT: 'pyramidingPolicy' defaults to 'mean'
- RETURNS: Exports a series of images for every n days from an as a GEE Asset
    which will show up as tasks in the 'Tasks' tab.
      * EXAMPLE: An image collection containing 52 images, one for each 
          week in a year, where you wish to return an image for every 4 weeks 
            * RESULT: Exports 13 images to the 'Tasks' tab with names of the form:
                assetId0, assetId4, assetId8, ..., assetId46, assetId50
*/


exports.exportImageCollectionToAsset = function(collection, assetId, roi, scale, crs, maxPixels, start, stepSize){
  print(collection);
  // the number of images in the image collection 
  var size = collection.size().getInfo();
  // return a list of the images in the collection
  var listOfImages = collection.toList(size);
  
  // collapse into a single, multi-band image
  var multiBandImg = ee.Image(listOfImages)
    .rename(ee.List.sequence(1, size).map(function (n) { return ee.String(n);}));

  // export each ith image as an Asset 
  Export.image.toAsset({
    image: multiBandImg,
    description: assetId,
    assetId: assetId,
    region: roi,
    scale: scale,
    crs: crs,
    maxPixels: maxPixels
  });
};
/*
  // iterate through the list of images starting at index 'start' by a step size of 'stepSiize'
  for (var i = start; i < size; i = i+stepSize) {
    // extract the ith image in the collection 
    var img = ee.Image(listOfImages.get(i));
    //var id = img.id().getInfo() || 'image_'+i.toString(); // line may not be necessary
    // create a name for each image exported: 'assetId' + 'image ID #'
    var id = ee.String(assetId).cat(ee.String('_fitted_day_')).cat(img.id()).getInfo();
    
    // export each ith image as an Asset 
    Export.image.toAsset({
      image: img,
      description: id,
      assetId: id,
      region: roi,
      scale: scale,
      crs: crs,
      maxPixels: maxPixels
    });
  }
};

*/

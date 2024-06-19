
/*
FUNCTION: exportImageCollectionToDrive
- Adapted from:
    * Author: Rodrigo E. Principe
    * License: Apache 2.0
    * email: fitoprincipe82@gmail.com
    https://gis.stackexchange.com/questions/407146/export-imagecollection-to-asset
- DESCRIPTION: Exports a series of images for every n days from an 
    ImageCollection as a TIF to Google Drive.
- PARAMETERS:
    * collection: (ImageCollection) an image collection
    * folder: (string) folder in Google Drive
    * assetId: (string) the name you choose for the exported Assets
    * rio: (Geometry.LinearRing|Geometry.Polygon|String) region of interest 
    * scale: (number) the scale or meters per pixel of the ImageCollection 
    * crs: (string) the CRS of the ImageCollection 
    * maxPixels: (number) Restrict the number of pixels in the export. By 
        default, you will see an error if the export exceeds 1e8 pixels
    * start: (number) starting index
    * stepSize: (number) step size
- PARAMETER DEFAULT: 'pyramidingPolicy' defaults to 'mean'
- RETURNS: Exports a series of images for every n days from an as a TIF
    which will show up as tasks in the 'Tasks' tab and be exported to Google
    Drive.
      * EXAMPLE: An image collection containing 30 images, one for each 
          day in a month, where you wish to return an image for every 4 days 
            * RESULT: Exports 7 images to the 'Tasks' tab with names of the form:
                assetIdDay0, assetIdDay4, assetIdDay8, ..., assetIdDay24, assetIdDay28
*/


exports.exportImageCollectionToDrive = function(collection, folder, assetId, roi, scale, crs, maxPixels, start, stepSize){
  // the number of images in the image collection 
  var size = collection.size().getInfo();
  // return a list of the images in the collection
  var listOfImages = collection.toList(size);

  // iterate through the list of images starting at index 0 by a step size of
    // size stepSiize
  for (var i = start; i < size; i = i+stepSize) {
    // extract the ith image in the collection 
    var img = ee.Image(listOfImages.get(i));
    //var id = img.id().getInfo() || 'image_'+i.toString(); // line may not be necessary
    // create a name for each image exported: 'assetId' + 'Day' + 'image ID #'
    var id = ee.String(assetId).cat('Day').cat(img.id()).getInfo();
    
  /// Export image to Drive ///
  /* NOTE: CANNOT HAVE PARAMETER VALUES IN BOTH 'scale' & 'dimensions'
      CAN ONLY HAVE ONE OR THE OTHER */
  Export.image.toDrive({
    folder: folder, 
    image: img, 
    description: id,
    scale: scale,
    maxPixels: maxPixels,
    crs: crs,
    region: roi 
    });
  }
};


// params and imports
var params = require('users/drewhart/seasonality/:params.js');

// get target projection CRS from SIF dataset
var SIF = ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN').first();
var SIF_proj = SIF.projection();
var SIF_scale = SIF_proj.nominalScale();

// calculate a global map, at our target resolution,
// of each pixel's minimum positive NIRv value
var calcNIRvPctDataAvail = function(nYears,
                                    endYear,
                                    dayStep
                                   ){
  // set dates for temporal filtering
  var end_yr = ee.Number(endYear);
  var start_yr = end_yr.subtract(nYears);
  var date_ending = ee.String('-01-01T00:00');
  var end_datestr = ee.String(end_yr).cat(date_ending);
  var start_datestr = ee.String(start_yr).cat(date_ending);
  var end_date = ee.Date(end_datestr);
  var start_date = ee.Date(start_datestr);
  // create list of target date strings, for slice-filtering by dayStep
  var total_day_span = end_date.difference(start_date,'day').round();
  var target_dates = ee.List.sequence(0, total_day_span, dayStep);
  var makeTargDateList = function(n) {
    // NOTE: format to match string format of system:index property
    return start_date.advance(n,'day').format('yyyy_MM_dd');
  };
  var targDateList = target_dates.map(makeTargDateList);
  // Load MODIS reflectance
  // NOTE: Band 2 contains NIR, 841-876 nm
  // NOTE: updated to MODIS v. 6.1, a la Reviewer 1's suggestion
  var NBARBands = ['Nadir_Reflectance_Band1', 'Nadir_Reflectance_Band2'];
  var MODISNBAR = ee.ImageCollection("MODIS/061/MCD43A4")
    // pull bands to be used for calculating NIRv
    .select(NBARBands)
    // immediately filter to only every dayStep-th image
    .filter(ee.Filter.inList('system:index', targDateList));
  var MODISNBARReproj = MODISNBAR
    .map(function(img){return img
      // and reproject to our target resolution
      .reduceResolution({reducer: ee.Reducer.mean(),
                         bestEffort: false,
                         maxPixels: 250});
    });
  // add NDVI and NIRv to that ImageCollection
  var NIRv = MODISNBARReproj.map(function(img){
      // calculate NDVI for each image
      // (NOTE: subtract 0.08 from NDVI, to attempt to
      // partially account for NDVI of bare soil per Badgley et al. 2017)
      // then calculate NIRv
      var NIRvImg = img.normalizedDifference(
        ['Nadir_Reflectance_Band2', 'Nadir_Reflectance_Band1'])
        .multiply(0)
        .add(1)
        .rename('NIRv')
        .copyProperties(img, ee.List(['system:time_start']));
      return NIRvImg;
    });
  var NIRvPctDataAvail = NIRv
    .map(function(img){return img.updateMask(img.gt(0))})
    // NOTE: FOR SOME REASON, DOING THIS USING ImageCollection.count() WAS TAKING
    //       FAR TOO LONG AND THEN FAILING BECAUSE OF OOM, EVEN FOR A SINGLE YEAR'S
    //       DATASET! NEEDED TO SWITCH TO SUM TO GET THIS TO RUN SUCCESSFULLY.
    .reduce(ee.Reducer.sum())
    .divide(MODISNBAR.size())
    .reproject({crs: SIF_proj, scale: SIF_scale});
  return NIRvPctDataAvail;
};

var NIRvPctDataAvail = calcNIRvPctDataAvail(params.NIRvNYears,
                                            params.NIRvEndYear,
                                            params.NIRvDayStep
                                           );

if (params.map_intermediates){
  print(NIRvPctDataAvail);
  print(NIRvPctDataAvail.projection().nominalScale().getInfo());
  Map.addLayer(NIRvPctDataAvail);
}

// export as a masking layer
Export.image.toAsset({image: NIRvPctDataAvail,
                      description: 'NIRv_ts_pct_data_availability',
                      assetId: 'users/drewhart/NIRv_ts_pct_data_availability',
                      region: params.roi,
                      pyramidingPolicy: 'mean',
                      maxPixels: params.maxPixels,
                      scale: SIF_scale.getInfo()
                    });

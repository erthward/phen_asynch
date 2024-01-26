//////////////////////
// ANCILLARY FUNCTIONS
//////////////////////


//------------
// format time
//------------

// use a universal time format that will work for both SIF and NIRv
// datasets, for both original and permuted data
// (function to be mapped to any dataset input to harmonic regression)
var reformatTimeStart = function(img){
  // NOTE: divide by 10^5, to cut down to a size that can be represented
  // without distortion by Int (no distortion because my datasets have a min
  // of 1-day temporal resolution, which is 24*60*60*1000 = 86,400,000 msec,
  // which will always remain a whole number after division by 10^5).
  // This way it can be fed into ImageCollection.remap as Int
  // (rather than Int64, which isn't valid),
  // then later I can multiply by 10^5 again to recover and use orig val
  var reformatted = ee.Number(ee.Image(img).get('system:time_start')).divide(100000);
  var out = ee.Image(img).set({'t': reformatted});
  return out;
};


//-----------
// mask water
//-----------

// NOTE: could be worth assessing for sensitivity to max occurrence of water
exports.maskWaterMODIS = function(imgColl, minWaterOccurr){
  // grab the JRC global surface water dataset; covers inland surface and coastal waters
  var waterMask1 = ee.Image("JRC/GSW1_0/GlobalSurfaceWater")
    .select('occurrence')
    .gt(ee.Image(minWaterOccurr))
    .unmask(0)
    .not();
  // grab a second MODIS water mask, which will cover more of the offshore waters
  var waterMask2 = ee.ImageCollection("MODIS/006/MOD44W")
  .select('water_mask')
  .reduce(ee.Reducer.sum())
  .gt(0)
  .not();
  // combine the two water masks
  var waterMask = waterMask1.and(waterMask2);
  // mask out water in all images in the collection
  var masked = imgColl.map(function(img){return img.updateMask(waterMask)});
  return masked;
};


/////////////////////////
// DATA-READING FUNCTIONS
/////////////////////////

//----------
// NIRv data
//----------

exports.getNIRvData = function(minWaterOccur, nYears){
  // Load MODIS reflectance (Band 2 should contain NIR, 841-876 nm)
  // NOTE: USING THESE DATES WILL ALWAYS END ON 12/31/18, THE LAST DATE OF THE SIF
  //       DATASET, BUT WILL EXTEND nYears BEFORE THAT
  var end_yr = ee.Number(2019);
  var start_yr = end_yr.subtract(nYears);
  var date_ending = ee.String('-01-01T00:00');
  var start_datestr = ee.String(start_yr).cat(date_ending);
  var end_datestr = ee.String(end_yr).cat(date_ending);
  var MODISReflectance = ee.ImageCollection("MODIS/006/MCD43A4")
    // first filter for a shorter date range, to hopefully help get around memory-usage issues
    // NOTE: for now taking a 10-year time slice to see how that does...
    .filterDate(start_datestr, end_datestr)
    // then reproject to the SIF resolution, using the mean reducer to aggregate,
    // to try to get around memory issues
    .map(function(img){return img.reduceResolution(ee.Reducer.mean(),false, 250)
    .reproject(ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN')
                    .first().projection())});
  // And add NDVI and NIRv to it
  var NIRv = MODISReflectance.map(function(img){
     var bands = ['Nadir_Reflectance_Band1', 'Nadir_Reflectance_Band2',
                   'nd', 'nd_1'];
      var new_bands = ['R', 'IR', 'NDVI', 'NIRv'];
      // add NDVI band to each image
      // NOTE: subtract 0.08 from NDVI, to attempt to partially account for NDVI of bare soil
      // (per Badgley et al. 2017)
      var NDVIAdded =  img
        .addBands(img.normalizedDifference(
        ['Nadir_Reflectance_Band2', 'Nadir_Reflectance_Band1']).subtract(0.08));
      // calculate NIRv (and mask out all values <=0, per Badgley et al. 2017
      var NIRvUnmasked = NDVIAdded.select('nd')
        .multiply(NDVIAdded.select('Nadir_Reflectance_Band2'));
      var ltezeroMask = NIRvUnmasked.gt(0);
      var NIRvAdded = NDVIAdded.addBands(NIRvUnmasked.mask(ltezeroMask));
      // select and rename the desired bands
      var NIRvImg = NIRvAdded
        .select(bands)
        .rename(new_bands);
      return NIRvImg;
    })
    // and reformat the start times
    .map(reformatTimeStart);
  // mask out pixels with water on greater than minWaterOccur pct of the time
  // in JRC surface-water dataset
  var NIRvmasked = exports.maskWaterMODIS(NIRv, minWaterOccur);
 
  return NIRvmasked;
};


//---------
// SIF data
//---------

// read in SIF data and reformat its time data
exports.getSIFData = function(minWaterOccur){
  var SIF = ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN')
    .map(reformatTimeStart);
  // mask out pixels with water on greater than minWaterOccur pct of the time
  // in JRC surface-water dataset
  var SIFmasked = exports.maskWaterMODIS(SIF, minWaterOccur);
  return SIFmasked;
};



//-------------
// climate data
//-------------

// load TerraClimate dataset (has multiple monthly climatological vars)
// and reformat its time data
exports.getTerraClimateData = function(minWaterOccur){
  
  // change "system:time_start" to the mean data for the monthly image
  // NOTE: NOT GREAT TO CREATE MIS-NAMED VARIABLE, BUT EASY WAY TO RECYCLE FUNCTION
  //       ABOVE AND QUICKLY PRODUCE A PROPERLY-NAMED VARIABLE
  var change_timestart_to_timemean = function(terraclim_img){
    var start = ee.Date(terraclim_img.get('system:time_start'));
    var end = ee.Date(terraclim_img.get('system:time_end'));
    var mean = ee.Date(terraclim_img.get('system:time_start'))
      .advance(end.difference(start, 'second').divide(2), 'second')
      // NOTE: convert back from Date to number, because 'system:time_start' must be a number
      .millis();
    return terraclim_img.set({'system:time_start': mean});
  };
  
  // load TerraClimate dataset (has multiple monthly climatological vars)
  // and reproject to our common SIF projection
  var SIF = ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN').first();
  var SIF_proj = SIF.projection();
  var SIF_scale = SIF_proj.nominalScale();      
  var terraclim = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
    .map(function(img){return img.reduceResolution({reducer: ee.Reducer.mean(),
                           bestEffort: false,
                           maxPixels: 250})
          .reproject({crs: SIF_proj, scale: SIF_scale})})
    .map(change_timestart_to_timemean)
    .map(reformatTimeStart)
    // and calculate tmmean as mean(tmmn, tmmx)
    .map(function(terraclim_img){return terraclim_img
      .addBands(terraclim_img.select('tmmx')
                  .add(terraclim_img.select('tmmn').rename('tmmx'))
                  .divide(2).rename('tmmean'), 
                ['tmmean']);
    });
  // mask out pixels with water on greater than minWaterOccur pct of the time
  // in JRC surface-water dataset
  var terraclim_masked = exports.maskWaterMODIS(terraclim, minWaterOccur);
  return terraclim_masked;
};
    
 
   
//-----------
// cloud data
//-----------

// load MODIS cloud data,
// and reformat its time data
exports.getMODISCloudData = function(minWaterOccur){
     
  // get start and end dates
  var start_date = ee.Date('2010-01-01');
  var end_date = ee.Date('2020-01-01');

  // alternate version of the reformatTimeStart function above
  // (because system:time_start is not a property in this data;
  //  instead, it only has system:index)
  var reformatTimeFromSystemIndex = function(img){
    var reformatted = ee.Date.parse('YYYY_MM_dd', img.get('system:index')).millis().divide(100000);
    var out = ee.Image(img).set('t',reformatted);
    return out;
  };
  
  // another alternate version, which slices the system:index string
  // after it was concatenated together with itself by the inner join
  var reformatTimeFromConcatSystemIndex = function(img){
    var reformatted = ee.Date.parse('YYYY_MM_dd', ee.String(img.get('system:index')).slice(0, 10))
                                    .millis().divide(100000);
    var out = ee.Image(img)
      .set('t', reformatted)
      .set('system:time_start', ee.Date(reformatted.multiply(100000)).millis());
    return out;
  };
  
  // QA-bit code adapted from https://spatialthoughts.com/2021/08/19/qa-bands-bitmasks-gee/
  // (itself adapted from https://gis.stackexchange.com/a/349401/5160
  var extractBit10 = function(input) {
    var fromBit = ee.Number(10);
    var toBit = ee.Number(10);
    var maskSize = ee.Number(1).add(toBit).subtract(fromBit);
    var mask = ee.Number(1).leftShift(maskSize).subtract(1);
    return input.rightShift(fromBit).bitwiseAnd(mask);
  };
  
  var load_modis_cloud = function(asset_id){ 
    var modis_cloud = ee.ImageCollection(asset_id)
      .filterDate(start_date, end_date)
      .map(function(img) {
        // only keep data with a number of observations per 1km pix is not 0
        // (since it appears from the User Guide that that flags bad data that is otherwise filled in with all 1s
        // in the QA bands)
        var img_masked = img.updateMask(img.select('num_observations_1km').gt(0));
        // use bit 10 to extract binary cloud state 
        var cloud = extractBit10(img_masked.select('state_1km')).eq(1)
          // rename proportion cloud band
          .rename('prop_cloud');
        return cloud;
      });
      return modis_cloud;
  };
  // load terra and aqua data
  var terra = load_modis_cloud("MODIS/006/MOD09GA");
  var terra = terra
    .map(function(img){return reformatTimeFromSystemIndex(img)});
  var aqua = load_modis_cloud("MODIS/006/MYD09GA");
  var aqua = aqua
    .map(function(img){return reformatTimeFromSystemIndex(img)});
  var filter = ee.Filter.equals({
  leftField: 't',
  rightField: 't'
  });

  // Get mean values for the inner-joined images (from the same dates)
  //var simpleJoin = ee.Join.saveAll({matchesKey: 'matchesKey', outer:true});
  var innerJoin = ee.Join.inner();
  var innerJoined = ee.ImageCollection(innerJoin.apply(terra, aqua, filter));
  var avg_cloud_prop = innerJoined.map(function(feature) {
    return ee.ImageCollection([feature.get('primary'), feature.get('secondary')])
      .reduce(ee.Reducer.mean())
      .rename('cloud_prop');
  });
  // add the t property again (but now have to take 10 chars of concatenated system:index prop)
  var avg_cloud_with_t = avg_cloud_prop
    .map(function(img){return reformatTimeFromConcatSystemIndex(img)});
  
  // last, reproject to SIF res
  var SIF = ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN').first();
  var SIF_proj = SIF.projection();
  var SIF_scale = SIF_proj.nominalScale();
  var avg_cloud_with_t_reproj = avg_cloud_with_t.map(function(img){
    return img.reproject({crs: SIF_proj, scale: SIF_scale})
      .reduceResolution({reducer: ee.Reducer.mean(),
                     bestEffort: false,
                     maxPixels: 250})
;
  });
  return avg_cloud_with_t_reproj;
};
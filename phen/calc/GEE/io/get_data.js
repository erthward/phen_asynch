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

exports.maskWaterMODIS = function(imgColl){
  // grab a second MODIS water mask, which will cover more of the offshore waters
  var waterMask = ee.ImageCollection("MODIS/006/MOD44W")
  .select('water_mask')
  .reduce(ee.Reducer.sum())
  .gt(0)
  .not();
  // mask out water in all images in the collection
  var masked = ee.ImageCollection(imgColl).map(function(img){return img.updateMask(waterMask)});
  return masked;
};


/////////////////////////
// DATA-READING FUNCTIONS
/////////////////////////

//----------
// NIRv data
//----------

exports.getNIRvData = function(maskLowBRDFAlbQual,
                               maxBRDFAlbQualVal,
                               maskWater,
                               nYears,
                               endYear,
                               dayStep,
                               clampToMinPos,
                               unmaskToMinPos,
                               maskLteZero
                              ){
  // get target projection CRS from SIF dataset
  var SIF = ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN').first();
  var SIF_proj = SIF.projection();
  var SIF_scale = SIF_proj.nominalScale();
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
  // load full-quality-info product and devise QA mask
  if (maskLowBRDFAlbQual){
    var QABands = ['BRDF_Albedo_Band_Quality_Band2', 'BRDF_Albedo_Band_Quality_Band1'];
    var nBands = QABands.length;
    var MCD43A2 = ee.ImageCollection("MODIS/061/MCD43A2")
      // get the bands needed for the QA mask
      .select(QABands)
      // immediately filter to only every dayStep-th image
      .filter(ee.Filter.inList('system:index', targDateList));
    // flag low quality pixels
    var QAMask = MCD43A2
      .map(function(img){return img.lte(maxBRDFAlbQualVal)
                                   .reduce(ee.Reducer.sum())
                                   .eq(nBands)});
  }
  // Load MODIS reflectance
  // NOTE: Band 2 contains NIR, 841-876 nm
  // NOTE: updated to MODIS v. 6.1, a la Reviewer 1's suggestion
  var NBARBands = ['Nadir_Reflectance_Band1', 'Nadir_Reflectance_Band2'];
  var MODISNBAR = ee.ImageCollection("MODIS/061/MCD43A4")
    // pull bands to be used for calculating NIRv
    .select(NBARBands)
    // immediately filter to only every dayStep-th image
    .filter(ee.Filter.inList('system:index', targDateList));
  if (maskLowBRDFAlbQual){
    var MODISNBARMasked = MODISNBAR
      // and mask with corresponding date from the QAMask ImageCollection
      // NOTE: the QAMask band is called 'sum' because of the output of the sum method
      .map(function(img){return img
        .updateMask(QAMask.filter(ee.Filter.eq('system:index', img.get('system:index'))).first().select('sum'))
      // and reproject to our target resolution
        .reduceResolution({reducer: ee.Reducer.mean(),
                           bestEffort: false,
                           maxPixels: 250})
        .reproject({crs: SIF_proj, scale: SIF_scale});
      });
  } else {
    var MODISNBARMasked = MODISNBAR
      .map(function(img){return img
        // and reproject to our target resolution
        .reduceResolution({reducer: ee.Reducer.mean(),
                           bestEffort: false,
                           maxPixels: 250})
        .reproject({crs: SIF_proj, scale: SIF_scale});
      });
  }
  // add NDVI and NIRv to that ImageCollection
  var NIRv = MODISNBARMasked.map(function(img){
      // calculate NDVI for each image
      // (NOTE: subtract 0.08 from NDVI, to attempt to
      // partially account for NDVI of bare soil per Badgley et al. 2017)
      // then calculate NIRv
      var NIRvImg = img.normalizedDifference(
        ['Nadir_Reflectance_Band2', 'Nadir_Reflectance_Band1'])
        .subtract(0.08)
        .multiply(img.select('Nadir_Reflectance_Band2'))
        .rename('NIRv')
        .copyProperties(img, ee.List(['system:time_start']));
      return NIRvImg;
    })
    // and reformat the start times by dividing by 10^5 and turning into a simpler 't' property
    .map(reformatTimeStart);
  // clamp zeros and negative values to their min observed NIRv values, if requested.
  // occurs predominantly over snowy regions and open water and barren areas;
  // snow, mostly in boreal winter, leads to pixel dropout when it persists long enough,
  // so this provides a reasonable backfill value that can retain those pixels and fit good LSP curves;
  // barren and open water get dropped during LC masking, so not a concern
  if (clampToMinPos | unmaskToMinPos){
    var minPosNIRv = ee.Image('users/drewhart/NIRv_min_pos');
  }
  if (clampToMinPos){
    var NIRv = NIRv.map(function(img){return img.where(img.lte(0), minPosNIRv)});
  } else if (maskLteZero) {
  // otherwise just drop all values <= 0
  // (unless maskLteZero if false, which is only the case for calculating monthly data-availability masks)
    var NIRv = NIRv.map(function(img){return img.updateMask(img.gt(0))});
  }
  // unmask all missing values to the min observed positive NIRv value, if requested
  // (just for exploring regression-fitting artefacts in places with long gaps)
  if (unmaskToMinPos){
    var NIRv = NIRv.unmask(minPosNIRv);
  }
  if (maskWater){
    // mask out water, if requested
    var NIRv = exports.maskWaterMODIS(NIRv);
  }
  return NIRv;
};


//---------
// SIF data
//---------

// read in SIF data and reformat its time data
exports.getSIFData = function(maskWater, clampToMinPos, unmaskToMinPos){
  var SIF = ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN')
    .map(reformatTimeStart);
  // clamp zeros and negative values to their min observed SIF values, if requested,
  // to match procedure applied to NIRv data (for straightforward comparison/assessment downstream)
  if (clampToMinPos | unmaskToMinPos){
    var minPosSIF = SIF.map(function(img){return img.updateMask(img.gt(0))}).reduce(ee.Reducer.min());
  }
  if (clampToMinPos){
    var SIF = SIF.map(function(img){return img.where(img.lte(0), minPosSIF)});
  }
  // otherwise don't worry about dropping all values <= 0, 
  // since we're not using SIF data except as a comparison to/assessment of NIRv results anyhow
  
  // unmask all missing values to the min observed positive SIF value, if requested,
  // to match procedure applied to NIRv data
  if (unmaskToMinPos){
    var SIF = SIF.unmask(minPosSIF);
  }
  // mask out water, if requested
  if (maskWater){
    var SIF = exports.maskWaterMODIS(SIF);
  }
  return SIF;
};



//-------------
// climate data
//-------------

// load TerraClimate dataset (has multiple monthly climatological vars)
// and reformat its time data
exports.getTerraClimateData = function(maskWater){
  
  // change "system:time_start" to the mean data for the monthly image
  // NOTE: NOT GREAT TO CREATE MIS-NAMED VARIABLE, BUT EASY WAY TO RECYCLE FUNCTION
  //       ABOVE AND QUICKLY PRODUCE A PROPERLY-NAMED VARIABLE
  var change_timestart_to_timemean = function(TerraClim_img){
    var start = ee.Date(TerraClim_img.get('system:time_start'));
    var end = ee.Date(TerraClim_img.get('system:time_end'));
    var mean = ee.Date(TerraClim_img.get('system:time_start'))
      .advance(end.difference(start, 'second').divide(2), 'second')
      // NOTE: convert back from Date to number, because 'system:time_start' must be a number
      .millis();
    return TerraClim_img.set({'system:time_start': mean});
  };
  
  // load TerraClimate dataset (has multiple monthly climatological vars)
  // and reproject to our common SIF projection
  var SIF = ee.ImageCollection('users/drewhart/seasonality_data/OCO2_SIF_ANN').first();
  var SIF_proj = SIF.projection();
  var SIF_scale = SIF_proj.nominalScale();      
  var TerraClim = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")
    .map(function(img){return img.reduceResolution({reducer: ee.Reducer.mean(),
                           bestEffort: false,
                           maxPixels: 250})
          .reproject({crs: SIF_proj, scale: SIF_scale})})
    .map(change_timestart_to_timemean)
    .map(reformatTimeStart)
    // and calculate tmmean as mean(tmmn, tmmx)
    .map(function(TerraClim_img){return TerraClim_img
      .addBands(TerraClim_img.select('tmmx')
                  .add(TerraClim_img.select('tmmn').rename('tmmx'))
                  .divide(2).rename('tmmean'), 
                ['tmmean']);
    });
  // mask out water, if requested
  if (maskWater){
    var TerraClim = exports.maskWaterMODIS(TerraClim);
  }
  return TerraClim;
};
    
 
   
//-----------
// cloud data
//-----------

// load MODIS cloud data,
// and reformat its time data
exports.getMODISCloudData = function(){
     
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
  var terra_imc = load_modis_cloud("MODIS/061/MOD09GA");
  var terra = terra_imc
    .map(function(img){return reformatTimeFromSystemIndex(img)});
  var aqua_imc = load_modis_cloud("MODIS/061/MYD09GA");
  var aqua = aqua_imc
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
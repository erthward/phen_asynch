// exports the raw NIRv data for a given set of sites
// (in a messy format, but that's no problem, very easy to clean up with Python,
// and I'm just flat out sick of working in GEE)

// load params used in main analysis
var params = require('users/drewhart/seasonality/:params.js');

// load data-prepping functions
var gDat = require('users/drewhart/seasonality/:io/get_data.js');                                                                                                                                                                                             

// load NIRv data
var NIRv = gDat.getNIRvData(params.maskLowBRDFAlbQual,
                            params.maxBRDFAlbQualVal,
                            params.maskWaterSeparately,
                            params.NIRvNYears,
                            params.NIRvEndYear,
                            params.NIRvDayStep,
                            params.NIRvClampToMinPos,
                            params.NIRvUnmaskToMinPos,
                            params.NIRvMaskLteZero
                           );
if (params.verbose){
  print('NIRv', NIRv);
}

// load the FLUXNET sites
var flux_sites = ee.FeatureCollection('users/drewhart/FLUXNET_sites');

// filter to the sites we want to plot
var site_names_to_extract = ['Blodgett Forest', 'Tonzi Ranch', 'Twitchell Wetland West Pond'];
var sites_to_extract = flux_sites.filter(ee.Filter.inList('SITE_NAME', site_names_to_extract));
if (params.verbose){
  print('sites_to_extract', sites_to_extract);
}

var pts_for_extraction = ee.FeatureCollection(sites_to_extract);    

var proj = ee.Image(NIRv.first()).projection();

var NIRv_pts = NIRv.map(function(img){
  var features = pts_for_extraction.map(function(f) {return f.set('date', img.get('date'))});
  return img.reduceRegions(features, ee.Reducer.mean(), 5, proj);
}).flatten();

if (params.verbose){
  print('NIRv_pts', NIRv_pts.limit(100));
}

// export result to Drive
Export.table.toDrive({'collection': NIRv_pts,
                      'description': 'raw_NIRv_at_CA_sites',
                      'folder': 'LSP_outputs_from_GEE',
                      'fileFormat': 'GeoJSON',
                      //'selectors': ['name', 'NIRv'],
});


// get parameters
var params = require('users/drewhart/seasonality/:params.js');
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
var SIF = gDat.getSIFData(false);
var proj = SIF.first().projection();
var scale = proj.nominalScale().getInfo();

// load and display Hansen 2022 LULCC dataset
var tiles = ee.ImageCollection('projects/glad/GLCmap2019');
var mosaic = tiles.mosaic().setDefaultProjection(tiles.first().projection());
var mosaicVizParam = {min:0,max:255,palette:['FEFECC','FDFDC8','FCFCC4','FAFAC0','F9F9BC','F7F7B8','F6F6B4','F4F4B0','F2F2AC','F1F1A8','EFEFA4',
'EEEEA1','ECEC9D','EBEB99','E9E995','E8E891','E6E68D','E5E589','E3E385','E2E281','E0E07D','DFDF7A','DDDD76','DBDB72','DADA6E','D8D86A','D5D564',
'D4D460','D2D25C','D1D158','CFCF54','CDCD50','CCCC4D','CACA49','C9C945','C7C741','C6C63D','C4C439','C3C335','C1C131','C0C02D','BEBE29','BDBD26',
'BBBB22','BABA1E','B8B81A','B6B616','B5B512','B3B30E','B2B20A','B0B006','609C60','5E9A5E','5C985C','5A965A','589458','569256','549054','528E52',
'508C50','4E8A4E','4C884C','4A864A','488448','468246','448044','427E42','407C40','3E7A3E','3C783C','3A763A','387438','367236','347034','326E32',
'316D31','2F6C2F','2C6A2C','296829','276627','246524','216321','1E611E','1B5E1B','175B17','145A14','115811','0E560E','0B540B','095309','065106',
'033303','FFABFF','FFA5FF','FF9EFF','FF98FF','FF91FF','FF8AFF','FF83FF','FF7DFF','FF76FF','FF6FFF','FF68FF','FF62FF','FF5AFF','FF53FF','FF4CFF',
'FF45FF','FF3EFF','FF38FF','FF31FF','FF2AFF','FF23FF','FF1CFF','FF16FF','FF0FFF','FF0000','000000','000000','000000','BFC0C0','BCBFC0','B8BEC2',
'B4BCC3','B1BBC4','ADBAC5','A9B9C6','A6B8C7','A2B7C9','9EB6CA','9AB5CB','97B4CC','93B2CD','90B2CE','8DB0CF','89AFD0','85AED1','82ADD3','7EACD4',
'7AABD5','77AAD6','73A9D7','70A8D8','6CA7D9','68A5DA','64A3DC','60A2DD','5CA0DE','589FDF','559EE0','519DE1','4E9CE2','4A9BE3','469AE4','4399E6',
'3F98E7','3B97E8','3895E9','3494EA','3193EB','2E92EC','2A91ED','2690EE','238FF0','1F8EF1','1C8DF2','188CF3','148BF4','118AF5','0D88F6','0986F7',
'9DC7C7','99C5C5','95C3C3','90C1C1','8CBFBF','87BDBD','83BBBB','7FB9B9','7AB7B7','76B5B5','72B3B3','6DB1B1','67AEAE','62ABAB','5EA9A9','5AA7A7',
'55A5A5','51A3A3','4DA1A1','489F9F','449D9D','3F9B9B','3A9999','369797','327C7C','327A7B','31787A','317678','307477','2F7275','2F7074','2F6E73',
'2C6A6F','2B696E','2B676D','2B656C','2A636A','296169','295F67','285D66','275B65','8C7BF0','8B77EF','8B74EF','8A71EE','8A6DEE','896AED','8967ED',
'8864EC','8861EC','875DEB','875AEB','8657EA','8351E7','834EE7','824BE6','8248E6','8144E5','8141E5','803EE4','803BE4','7F38E3','7F34E3','7E31E2',
'7D2EE1','C80000','000000','000000','000000','00F4F4','00E8E8','00DDDD','00D0D0','00C5C5','00B7B7','00ACAC','009F9F','009494','008888','00007D',
'FFFFFF','FF7D00','000000','000032','010101']};
if (params.map_intermediates){
  Map.addLayer(mosaic, mosaicVizParam,'global land cover and land use', false);
}

// reclass LULCC map
var lulcc = mosaic
  .gte(92).and(mosaic.lte(116))               // terra firma tree gain or loss 
  .or(mosaic.gte(212).and(mosaic.lte(236)))   // wetland tree gain or loss
  .or(mosaic.gte(240).and(mosaic.lte(249)))   // built-up land
                                              // (NOTE: these pixel values could actually be multiplied by the fractions
                                              //        in their bin definitions, but built-up land has a small footprint
                                              //        should already be largely excluded from our asynchrony map,
                                              //        so this is unlikely to make any real difference in the rf model
  .or(mosaic.eq(252));                        // cropland
if (params.map_intermediates){
  Map.addLayer(lulcc, {min:0, max:1, palette:['white', 'red']}, 'land use and land cover change', true);
}

// aggregate percent subpixel lulcc to 1km and export
var lulccPct_1kagg = lulcc
  .reduceResolution({reducer: ee.Reducer.mean(),
                     maxPixels: 2000})
  .reproject({crs: proj, scale: 1000})
  .rename('pct_subpix_lulcc');
Export.image.toAsset({
  image: lulccPct_1kagg,
  description: 'lulccPct_1kagg',
  assetId: 'users/drewhart/lulccPct_1kagg',
  region: params.roi,
  pyramidingPolicy: 'mean',
  maxPixels: params.maxPixels,
  scale: 1000,
});

// read in that 1km-res aggregate and aggregate again to our target 5km res
var lulccPct_1kagg = ee.Image('users/drewhart/lulccPct_1kagg');
if (params.map_intermediates){
  Map.addLayer(lulccPct_1kagg, {min:0, max:1, palette:['white', 'red']}, 'pct subpixel lulcc 1kagg', true);
}
var lulccPct_5kagg = lulccPct_1kagg
  .reduceResolution({reducer: ee.Reducer.mean(),
                     maxPixels: 2000})
  .reproject({crs: proj, scale: scale})                   
  .rename('pct_subpix_lulcc');
if (params.map_intermediates){
  Map.addLayer(lulccPct_5kagg, {min:0, max:1, palette:['white', 'red']}, 'pct subpixel lulcc 5kagg', true);
}

// mask to only pixels whose LSP signal is being factored into the LSP asynchrony map
var strict_mask = ee.Image('users/drewhart/LSP_mask_STRICT');
var lulccPct_mask = lulccPct_5kagg.mask(strict_mask);
if (params.map_intermediates){
  Map.addLayer(lulccPct_mask, {palette: ['white', 'red'], min:0, max:1}, 'pct subpixel lulcc 5kagg, strict-masked'); 
}

// calculate neighborhood mean fraction of LULCC within 100 km neighborhood
var lulccPct_neighMean = lulccPct_mask
  .reduceNeighborhood({reducer: ee.Reducer.mean(),
                       kernel: ee.Kernel.circle(100000, 'meters'),
                      });
Map.addLayer(lulccPct_neighMean, {min:0, max:1, palette:['white', 'red']}, 'neigh mean LULCC fraction');

// export both the 5km agg map and its 100 km-neighborhood mean to file
Export.image.toDrive({
  image: lulccPct_mask,
  description: 'hansen_lulcc_pct',
  folder: 'LSP_asynch_drivers_analysis_covars',
  scale: scale,
  region: params.roi,
  fileFormat: 'GeoTIFF',
  maxPixels: params.maxPixels,
  shardSize: params.shardSize,
  fileDimensions: params.fileDimensions,
});
Export.image.toDrive({
  image: lulccPct_neighMean,
  description: 'hansen_lulcc_pct_neigh_mean',
  folder: 'LSP_asynch_drivers_analysis_covars',
  scale: scale,
  region: params.roi,
  fileFormat: 'GeoTIFF',
  maxPixels: params.maxPixels,
  shardSize: params.shardSize,
  fileDimensions: params.fileDimensions,
});



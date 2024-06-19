var palettes = require('users/gena/packages:palettes');

var esa_lc = ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V-C3/Global");
Map.addLayer(esa_lc.first().select('discrete_classification'),
             {palette: palettes.colorbrewer.Paired[12], opacity:0.6}, 'ESA LC', false);

var ecoregions = ee.FeatureCollection("RESOLVE/ECOREGIONS/2017");
Map.addLayer(ecoregions, {opacity:0.3}, 'ecoregion boundaries', false);
Map.addLayer(ecoregions.reduceToImage(['BIOME_NUM'], ee.Reducer.first()),
             {palette: palettes.colorbrewer.Set3[12], 
              min: 0, max: 14, opacity:0.6},
             'biomes', false);

var eofs = ee.Image('users/drewhart/NIRv_phen_EOFs');
// feed EOF2 into red and EOF1 into blue (leave green empty)
var eofs = eofs.select(['b1'])
  .addBands(eofs.select(['b2']))
  .addBands(eofs.select(['b3']));
print('EOFs', eofs);
Map.addLayer(eofs, {opacity:0.75}, 'EOFs');

var asynch = ee.Image('users/drewhart/NIRv_global_asynch');
Map.addLayer(asynch.select('b3'),
             {min:-0.00002, max:0.0002,
               palette: palettes.matplotlib.magma[7], opacity:0.75}, 'asynchrony', false);
               
Map.setOptions('SATELLITE');


//-------------------------
// params for model fitting
//-------------------------

var n = 2;
var harmonics = ee.List.sequence(1, n);


//--------
// imports
//--------

var hReg = require('users/drewhart/seasonality/:fft/harmonic_regression.js');
var gDat = require('users/drewhart/seasonality/:io/get_data.js');
var cbar = require('users/drewhart/seasonality/:viz/make_colorbar.js');

//---------
// palettes
//---------

// load palettes
var palettes = require('users/gena/packages:palettes'); 


//----------
// NIRv data
//----------

var NIRv = gDat.getNIRvData(50, 10);
//print('NIRv', NIRv.first());


//---------
// SIF data
//---------

// read in SIF data and reformat its time data
var SIF = gDat.getSIFData(50);
//print('SIF', SIF);



//----------------------
// run regressions: NIRv
//----------------------

// calculate regression with 4 freqs
var dependent = 'NIRv';
var bandsForViz = ee.List(['sin_1', 'cos_1']);
var reg = ee.ImageCollection(hReg.calcHarmonicRegression(NIRv, dependent,
                                                         harmonics, bandsForViz));
print('regression', reg.first());
var ts = ee.ImageCollection(reg.toList(365*4,1)).select(['NIRv', 'fitted']);
print('fitted time series', ts);

//-------
// map it
//-------
Map.addLayer(reg.first().select('phase'),
             {palette: palettes.kovesi.cyclic_mygbm_30_95_c78_s25[7], 
              min:0, max:1, opacity: 0.8}, 'phase', false);
//Map.add(cbar.makeLegend('annual-frequency phase', 'Jan 1', 'Jul 1', 'Dec 31',
//        palettes.kovesi.cyclic_mygbm_30_95_c78_s25[7], 'top-left'));



//------------------------
// point comparison charts
//------------------------

var getTimeSeries = function(img_coll, pt){
  var timeSeriesColl = img_coll.map(function(img){
    var reduced = img.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: pt,
      scale: 5000 // TODO: How can I figure out the native res of the SIF asset?
                 //       And does it matter at all if I stipulate a scale less than the native res?
    }).get('fitted');
    // TODO: REPLACE THE APPROACH IN THE NEXT LINE WITH EITHER SOME JUSTIFIABLE INTERPOLATION
    // OR A MUTUAL FILTERING OF MISSING VALUES
    var reduced = ee.Algorithms.If(ee.Algorithms.IsEqual(reduced, null), 0, reduced);
    return img.set('extracted', reduced);
  }, false);
  var timeSeries = timeSeriesColl.aggregate_array('extracted');
  return timeSeries;
};


//------------------
// POINT-DRAWING APP
//------------------

var drawingTools = Map.drawingTools();
drawingTools.setShown(false);
while (drawingTools.layers().length() > 0) {
  var layer = drawingTools.layers().get(0);
  drawingTools.layers().remove(layer);
}
var dummyGeometry =
    ui.Map.GeometryLayer({geometries: null, name: 'geometry', color: '23cba7'});
drawingTools.layers().add(dummyGeometry);

function clearGeometries() {
  while (drawingTools.layers().length() > 0) {
    var layer = drawingTools.layers().get(0);
    drawingTools.layers().remove(layer);
  }
}

function drawPoint() {
  drawingTools.setShape('point');
  drawingTools.draw();
}

//var scatchartPanel = ui.Panel({
//  style:
//      {height: '235px', width: '600px', position: 'top-right', shown: false}
//});
//Map.add(scatchartPanel);
var p1chartPanel = ui.Panel({
  style:
      {height: '235px', width: '600px', position: 'bottom-left', shown: false}
});
Map.add(p1chartPanel);
var p2chartPanel = ui.Panel({
  style:
      {height: '235px', width: '600px', position: 'bottom-right', shown: false}
});
Map.add(p2chartPanel);


function chartTimeSeries() {
  // Make the chart panel visible the first time a geometry is drawn.
  //if (!scatchartPanel.style().get('shown')) {
  //  scatchartPanel.style().set('shown', true);
  //}
   if (!p1chartPanel.style().get('shown')) {
    p1chartPanel.style().set('shown', true);
  }
   if (!p2chartPanel.style().get('shown')) {
    p2chartPanel.style().set('shown', true);
  }

  // Get the drawn geometry; it will define the reduction region.
  var pts = drawingTools.layers().get(0).getEeObject();
  var p1 = ee.Geometry.Point(pts.coordinates().get(0));
  var p2 = ee.Geometry.Point(pts.coordinates().get(1));
  
  // Get the time series at both points
  var ts1 = getTimeSeries(ts.select('fitted'), p1);
  var ts2 = getTimeSeries(ts.select('fitted'), p2); 

  // Set the drawing mode back to null; turns drawing off.
  drawingTools.setShape(null);

  // Chart time series for the selected area of interest.
  var p1chart = ui.Chart.image.series(ts.select(['NIRv', 'fitted']), p1)
    .setOptions({title: 'p1',
                 hAxis: {'title': ee.String(p1.coordinates())}});
  var p2chart = ui.Chart.image.series(ts.select(['NIRv', 'fitted']), p2)
    .setOptions({title: 'p2',
                 hAxis: {'title': p2.coordinates()}});
  //var scatchart = ui.Chart.array.values(ts2, 0, ts1)
  //  .setSeriesNames(['p2'])
  //  .setOptions({
  //    title: 'correlation between points 1 and 2',
  //    hAxis: {'title': 'p1'},
  //    vAxis: {'title': 'p2'},
  //    pointSize: 3,
//});


  // Replace the existing chart in the chart panel with the new chart.
  p1chartPanel.widgets().reset([p1chart]);
  p2chartPanel.widgets().reset([p2chart]);
  //scatchartPanel.widgets().reset([scatchart]);
}

drawingTools.onDraw(ui.util.debounce(chartTimeSeries, 500));
drawingTools.onEdit(ui.util.debounce(chartTimeSeries, 500));

var symbol = {
  pointDraw: '✍',
};

var controlPanel = ui.Panel({
  widgets: [
    ui.Label('Draw 2 points'),
    ui.Button({
      label: symbol.pointDraw + ' Draw a Point',
      onClick: drawPoint,
      style: {stretch: 'horizontal'}
    }),
    ui.Label('Wait for charts to render.'),
    ui.Button({
      label: 'Clear All Points',
      onClick: clearGeometries,
      style: {stretch: 'horizontal'}
    }),
  ],
  style: {position: 'top-left'},
  layout: null,
});
Map.add(controlPanel);


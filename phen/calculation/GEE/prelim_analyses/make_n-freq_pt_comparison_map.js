/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var geometry = /* color: #23cba7 */ee.Geometry.Point([-121.98961353746053, 40.484629182643616]);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
//-------
// params
//-------

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

var NIRvP = gDat.getNIRvData(50);
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
var dependent = 'b1';
var bandsForViz = ee.List(['sin_1', 'cos_1']);
var reg = ee.ImageCollection(hReg.calcHarmonicRegression(SIF.limit(1000), dependent,
                                                         harmonics, bandsForViz));
print('regression', reg.first());
var ts = ee.ImageCollection(reg.toList(101, 1)).select(['b1', 'fitted']);
print('fitted time series', ts);

//-------
// map it
//-------
Map.setOptions('TERRAIN');
Map.addLayer(reg.first().select('phase'),
             {palette: palettes.kovesi.cyclic_mygbm_30_95_c78_s25[7], 
              min:0, max:1, opacity: 0.8}, 'phase');
Map.add(cbar.makeLegend('annual-frequency phase', 'Jan 1', 'Jul 1', 'Dec 31',
        palettes.kovesi.cyclic_mygbm_30_95_c78_s25[7], 'top-left'));



//------------------------
// point comparison charts
//------------------------

var getTimeSeries = function(img_coll, pt){
  var timeSeriesColl = img_coll.map(function(img){
    var reduced = img.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: pt,
      scale: 500 // TODO: How can I figure out the native res of the SIF asset?
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

/*
var p1 = ee.Geometry.Point(-59.598, 6.979);
var p2 = ee.Geometry.Point(-59.949, 4.049);
print('p1', p1.coordinates());
print('p2', p2.coordinates());
Map.addLayer(p1, {color: 'FF0000'}, 'p1');
Map.addLayer(p2, {color: 'FF0000'}, 'p2');

var ts1 = getTimeSeries(ts.select('fitted'), p1);
var ts2 = getTimeSeries(ts.select('fitted'), p2); 

// plot sampled features as two time series and a scatter chart

print(ui.Chart.image.series(ts.select(['NIRvP', 'fitted']), p1)
  .setOptions({title: 'p1'}));

print(ui.Chart.image.series(ts.select(['NIRvP', 'fitted']), p2)
  .setOptions({title: 'p2'}));

var chart = ui.Chart.array.values(ts2, 0, ts1)
    .setSeriesNames(['p2'])
    .setOptions({
      title: 'correlation between points 1 and 2',
      hAxis: {'title': 'p1'},
      vAxis: {'title': 'p2'},
      pointSize: 3,
});
print(chart);
*/





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

var scatchartPanel = ui.Panel({
  style:
      {height: '235px', width: '600px', position: 'top-right', shown: false}
});
Map.add(scatchartPanel);
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


function chartSIFTimeSeries() {
  // Make the chart panel visible the first time a geometry is drawn.
  if (!scatchartPanel.style().get('shown')) {
    scatchartPanel.style().set('shown', true);
  }
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

  // Chart SIF time series for the selected area of interest.
  var p1chart = ui.Chart.image.series(ts.select(['b1', 'fitted']), p1)
    .setOptions({title: 'p1',
                 hAxis: {'title': ee.String(p1.coordinates())}});
  var p2chart = ui.Chart.image.series(ts.select(['b1', 'fitted']), p2)
    .setOptions({title: 'p2',
                 hAxis: {'title': p2.coordinates()}});
  var scatchart = ui.Chart.array.values(ts2, 0, ts1)
    .setSeriesNames(['p2'])
    .setOptions({
      title: 'correlation between points 1 and 2',
      hAxis: {'title': 'p1'},
      vAxis: {'title': 'p2'},
      pointSize: 3,
});


  // Replace the existing chart in the chart panel with the new chart.
  p1chartPanel.widgets().reset([p1chart]);
  p2chartPanel.widgets().reset([p2chart]);
  scatchartPanel.widgets().reset([scatchart]);
}

drawingTools.onDraw(ui.util.debounce(chartSIFTimeSeries, 500));
drawingTools.onEdit(ui.util.debounce(chartSIFTimeSeries, 500));

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



////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   STUFF FOR LAB MEETING   /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

// add 50km GEE-calculated asynch and year-fitted seasonality, to be able to inspect and compare


var pi = ee.Number(3.141592653589793);
var msecPerDay = ee.Number(86400000);


var makeSingleYearHarmonicImgColl = function(datasetProj){
  // get the integer days, 0 to 364
  var days = ee.List.sequence({start:0, end:364, step:1});
  // re-express those as radians in both the annual and semiannual freqs
  var annualRads = days.map(function(n){
    return ee.Number(n).divide(365).multiply(pi.multiply(2))});
  var semianRads = days.map(function(n){
    return ee.Number(n).divide(365).multiply(pi.multiply(4))});
  // create ImageCollections of constant images of the sines and cosines
  // of those radian-transformed values
  var annualSin = ee.ImageCollection(annualRads.map(function(n){
    return ee.Image.constant(n).sin().rename('ann_sin')}));
  var annualCos = ee.ImageCollection(annualRads.map(function(n){
    return ee.Image.constant(n).cos().rename('ann_cos')}));
  var semianSin = ee.ImageCollection(semianRads.map(function(n){
    return ee.Image.constant(n).sin().rename('sem_sin')}));
  var semianCos = ee.ImageCollection(semianRads.map(function(n){
    return ee.Image.constant(n).cos().rename('sem_cos')}));
  // combine all into one ImageCollection
  var year = ee.ImageCollection(ee.ImageCollection(ee.ImageCollection(annualSin
    .combine(annualCos))
    .combine(semianSin))
    .combine(semianCos))
  // reproject to the datasetProj
    .map(function(img){return img.reproject(datasetProj)});
  var yearWithIndex = year.map(function(img){
    return img.set({'num_index': ee.Number.parse(ee.Image(img).get('system:index')),
      // 'system:time_start': ee.Date(ee.Number.parse(ee.Image(img).get('system:index')).multiply(msecPerDay))})})
    })})
    .sort('num_index');
  return yearWithIndex;
};




// function that takes the first image from a harmonic regression ImgColl
// (the image which contains all the coefficients and so on)
// and a year's worth of annual and semiannual-freq sine and cosine harmonics,
// and returns a new 365-length ImgColl of the fitted values for each day of the year
var fitSingleYearDailyValues = function(harmRegResults, proj){
  var year = makeSingleYearHarmonicImgColl(proj);
  var fittedYear = year.map(function(img){
    return harmRegResults.select('constant')
      .add(harmRegResults.select('sin_1').multiply(img.select('ann_sin')))
      .add(harmRegResults.select('cos_1').multiply(img.select('ann_cos')))
      .add(harmRegResults.select('sin_2').multiply(img.select('sem_sin')))
      .add(harmRegResults.select('cos_2').multiply(img.select('sem_cos')))
      .rename('fitted')
      // copy the properties, to carry over the num_index property I created
      .copyProperties(img)});
  return fittedYear;
};


var myPal = ['#000000', '#171616', '#362f2f', '#523a3a', 
               '#7d3e3e', '#9e3737', '#c42727', '#f52525'];
var asynchrony50 = ee.Image("users/drewhart/masked_asynch_SIF_TEST_50000m_neighhood");
Map.addLayer(asynchrony50, {min:4.4584500616377436e-8,
                            max: 0.000025863525487424156,
                            opacity:0.8, palette: myPal }, 'asynchrony 50km');

var fittedYear = fitSingleYearDailyValues(reg.first(), asynchrony50.projection());
Map.addLayer(fittedYear, {opacity:0.0001}, 'fittedYear');


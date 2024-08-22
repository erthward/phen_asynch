// This script creates a GEE app that provides access to the main mapping
// products from Terasaki Hart et al. 2024. (URL for the paper:
// DOI for the paper). This work maps global estimates of the diversity
// and asynchrony of annual vegetation seasonality
// (also known as 'land surface phenology' (LSP))
// and explores its ecological drivers and evolutionary implications.
//
//
// TO USE:
//   1. Click the 'Run' button, in bar above this panel.
//   2. Once the map has loaded, use the 'Layers' button, in the map panel,
//      to add or remove map layers, or to control their transparency, colors, etc.
//   3. Once the map has loaded, click the '🖍 Draw points' button, then
//      point and click to drop markers at any two locations.
//   4. Wait for the line plots to render, providing a visualization of their raw
//      modeled annual land surface phenology patterns
//      (based on MODIS Near Infrared Reflectance of Vegetation; NIRv).
//   5. Click on markers, then click and drag to new locations, to dynamically
//      update line plots.
//   6. Click the '🗑 ️Clear all points' button, if needed,
//      then repeat the above process to create new points.
//
//
// All data are free to use and redistribute, with proper
// attribution under the MIT license and with proper
// citation of the original scientific paper (DOI ).
//
//=======================================================================================
//=======================================================================================
//=======================================================================================

////////////////////////////////////////////////////////////////////
// NIRv DATA-READING FUNCTIONS
////////////////////////////////////////////////////////////////////

// use a universal time format that will work for both SIF and NIRv
// datasets, for both original and permuted data
// (function to be mapped to any dataset input to harmonic regression)
var reformatTimeStartForViz = function(img){
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


// read the 20-year NIRv dataset
// (NOTE: this version is created to be used to load the NIRv data
//        for the data viewer, with values hard-coded to the same values used
//        for the analysis in the paper, thus reducing dependency on importing
//        of additional scripts)
var getNIRvDataForViz = function(){
  // define vars that are input args in the analysis-ready form of the function
  var endYear = 2021;
  var nYears = 20;
  var dayStep = 4;
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
    return start_date.advance(n, 'day').format('yyyy_MM_dd');
  };
  var targDateList = target_dates.map(makeTargDateList);
  // Load MODIS reflectance
  // NOTE: Band 2 contains NIR, 841-876 nm
  // NOTE: updated to MODIS v. 6.1, a la Reviewer 1's suggestion
  var NBARBands = ['Nadir_Reflectance_Band1', 'Nadir_Reflectance_Band2'];
  var MODISNBAR = ee.ImageCollection("MODIS/061/MCD43A4")
    // pull bands to be used for calculating NIRv
    .select(NBARBands)
    // immediately filter to only every dayStepth image
    .filter(ee.Filter.inList('system:index', targDateList));
  var MODISNBARReproj = MODISNBAR
    .map(function(img){return img
      // and reproject to our target resolution
      .reduceResolution({reducer: ee.Reducer.mean(),
                         bestEffort: false,
                         maxPixels: 250})
      .reproject({crs: SIF_proj, scale: SIF_scale});
    });
  // add NDVI and NIRv to that ImageCollection
  // (clamping NIRv to minimum positive values observed across time series)
  // (and properly rescaling all values by dividing by 10k, since we are not outputting
  // files we need to reduce the size of, and the displayed reflectance values must be correct)
  var minPosNIRv = ee.Image('users/drewhart/NIRv_min_pos').divide(10000);
  var NIRv = MODISNBARReproj.map(function(img){
      // calculate NDVI for each image
      // (NOTE: subtract 0.08 from NDVI, to attempt to
      // partially account for NDVI of bare soil per Badgley et al. 2017)
      // then calculate NIRv
      var NIRvImg = img.normalizedDifference(
        ['Nadir_Reflectance_Band2', 'Nadir_Reflectance_Band1'])
        .subtract(0.08)
        .multiply(img.select('Nadir_Reflectance_Band2'))
        .divide(10000)
        .rename('NIRv')
        .copyProperties(img, ee.List(['system:time_start']));
      return NIRvImg;
    })
    // and reformat the start times by dividing by 10^5 and turning into a simpler 't' property
    .map(reformatTimeStartForViz);
  var NIRv = NIRv.map(function(img){return img.where(img.lte(0), minPosNIRv)});
  return NIRv;
};




////////////////////////////////////////////////////////////////////
// HARMONIC REGRESSION FUNCTIONS
////////////////////////////////////////////////////////////////////

// Function to get a sequence of band names for harmonic terms.
var constructBandNames = function(base, list) {
  return ee.List(list).map(function(i) {
    return ee.String(base).cat(ee.Number(i).int());
  });
};


// Function to add a time band.
var addDependents = function(img) {
  // Compute time in fractional years since the epoch.
  var unreformattedTime = ee.Date(ee.Number(ee.Image(img).get('t')).multiply(100000));
  var years = unreformattedTime.difference('2000-01-01', 'year');
  var timeRadians = ee.Image(years.multiply(2 * Math.PI)).rename('t');
  var constant = ee.Image(1);
  return ee.Image(img).addBands(constant).addBands(timeRadians.float());
};

// Function to compute the specified number of harmonics
// and add them as bands.  Assumes the time band is present.
var addHarmonics = function(freqs, sinNames, cosNames) {
  return function(img) {
    // Make an image of frequencies.
    var frequencies = ee.Image.constant(freqs);
    // This band should represent time in radians.
    var time = ee.Image(img).select('t');
    // Get the cosine terms.
    var cosines = time.multiply(frequencies).cos().rename(cosNames);
    // Get the sin terms.
    var sines = time.multiply(frequencies).sin().rename(sinNames);
    return ee.Image(img).addBands(cosines).addBands(sines);
  };
};

// function to compute R2s and return as an Image
// (assumes an ImageCollection with Images each having a band of fitted values 
// with the name 'fitted' and a band of original dependent-variable values with
// the name being whatever value the argument `dependent` holds)
var calcR2 = function(imgColl, dependent){
  return imgColl
    // mask out fitted values on dates with missing input values
    // (so that R2 values should not exceed 1)
    .map(function(img){return img.updateMask(img.select(dependent).mask())})
    .select('fitted')
    .reduce(ee.Reducer.variance())
    .divide(imgColl
      .select(dependent)
      .reduce(ee.Reducer.variance()))
    .rename('R2');
};


// a version of the harmonic regression function, trimmed down solely for use for the data viewer
var calcHarmonicRegressionForViz = function(NIRvImgColl){
  
  // define vars that are input args in the analysis-ready form of the function
  var dependent = 'NIRv';
  var harmonics = ee.List.sequence(1, 2);
  
  //------------------------------------
  // Create the harmonic ImageCollection
  //------------------------------------
  
  // Construct lists of names for the harmonic terms.
  var cosNames = constructBandNames('cos_', harmonics);
  var sinNames = constructBandNames('sin_', harmonics);
  // Independent variables.
  var independents = ee.List(['constant', 't'])
    .cat(cosNames).cat(sinNames);
    
  // Add variables
  var harmonicImgColl = ee.ImageCollection(NIRvImgColl
    // add the dependent and harmonic variables
    .map(addDependents)
    .map(addHarmonics(harmonics, sinNames, cosNames))
    // and add independents as a property
    .set({'independents': independents}));
  
  //-------------------------
  // Calculate the regression
  //-------------------------
  
  // Calculate regression reduction and get output (a 4x1 array image)
  var harmonicTrend = harmonicImgColl
    .select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));
  
  // Turn the array image into a multi-band image of coefficients.
  var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);
  
  // calculate the residuals
  var harmonicTrendResiduals = harmonicTrend.select('residuals')
    .arrayProject([0])
    .arrayFlatten([['residuals']]);
    
  // Pull out the bands to be visualized
  var sin = harmonicTrendCoefficients.select(ee.String('sin_1'));
  var cos = harmonicTrendCoefficients.select(ee.String('cos_1'));
  
  // Do some math to turn a first-order Fourier model into
  // hue, saturation, and value in the range[0,1].
  var amplitude = cos.hypot(sin).multiply(5).rename('amplitude');
  var phase = sin.atan2(cos).unitScale(-Math.PI, Math.PI).rename('phase');
  var val = harmonicImgColl.select(dependent).reduce('mean');
  
  // Turn the HSV data into an RGB image and add it to the map.
  var seasonality = ee.Image.cat(phase, amplitude, val).hsvToRgb()
                            .rename(['seasonality_R',
                                     'seasonality_G',
                                     'seasonality_B']);
    
  // get SD of seasonality phase within sliding-window neighborhoods
  var phaseSD = phase.reduceNeighborhood({
    reducer: ee.Reducer.stdDev(),
    kernel: ee.Kernel.circle({radius: 100000, units: 'meters'})
  }).rename('phaseSD');
                           
  //---------------------------------
  // Calculate the fitted vals and R2
  //---------------------------------
  
  var fittedHarmonic = harmonicImgColl.map(function(img) {
    return ee.Image(img).addBands(
      img.select(independents)
         .multiply(harmonicTrendCoefficients)
         .reduce('sum')
         .rename('fitted'));
  });

  // Fix the time (in radians) of the maximum fitted value
  var tmax = fittedHarmonic.select(['fitted', 't'])
    .reduce(ee.Reducer.max(2))
    .rename('max_fitted', 't_max')
    .select('t_max')
    .mod(6.283185307179586);
    
  // Calculate the R2s
  var R2 = calcR2(fittedHarmonic, dependent);
  
  // compose final output ImageCollection for the model
  var summary = seasonality.addBands(harmonicTrendCoefficients)
                           .addBands(harmonicTrendResiduals)
                           .addBands(tmax)
                           .addBands(R2)
                           // add max possible time_start val as property, to make it
                           // easy later to separate from original NIRvImgColl it is merged with
                           .set({'system:time_start': 999999999999});
  
  // merge the summary image and the rest of the ImageCollection
  var harmonicReg = ee.ImageCollection(summary)
    .merge(fittedHarmonic);
  // extract just the raw and fitted NIRv values
  var output = ee.ImageCollection(harmonicReg.toList(NIRvImgColl.size(),1)).select(['NIRv', 'fitted']);
  return output;
};



/////////////////////////////////////////////////////////
// UI FUNCTIONS
/////////////////////////////////////////////////////////

var palettes = require('users/gena/packages:palettes');

var ColorBar = function(palette) {
  return ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: {
      bbox: [0, 0, 1, 0.1],
      dimensions: '200x20',
      format: 'png',
      min: 0,
      max: 1,
      palette: palette,
    },
    style: {stretch: 'horizontal', margin: '0px 8px'},
  });
};

var makeLegend = function(title, low, mid, high, palette, pos) {
  var titlePanel = ui.Panel(
      [
        ui.Label(title, {textAlign: 'center', stretch: 'horizontal'})
      ],
      ui.Panel.Layout.flow('horizontal'));
  var labelPanel = ui.Panel(
      [
        ui.Label(low, {margin: '4px 8px'}),
        ui.Label(mid, {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
        ui.Label(high, {margin: '4px 8px'})
      ],
      ui.Panel.Layout.flow('horizontal'));
  var panel = ui.Panel([titlePanel, ColorBar(palette), labelPanel]);
  panel.style().set({
      width: '300px',
      position: 'top-left'
    });
  return panel;
};

// get palette to use for asynchrony
// (because calling the 'reverse' method more than once reverses the palette again!)
var asynchPalette = palettes.crameri.lajolla[25].reverse();

// create an asynchrony legend (to be added when that layer is displayed, otherwise removed)
var asynchLegend = makeLegend('LSP asynchrony',
                              'low',
                              'mod',
                              'high',
                              asynchPalette
                             );

// function to set up and stylize map
var createMap = function(){
  // clear existing map
  ui.root.clear();
  // create the map object and add it
  var map = ui.Map();
  ui.root.add(map);
  // Set basemap options.
  map.setOptions('HYBRID');

  // Set visibility options to remove geometry creator, map type controller,
  // and layer list.
  map.setControlVisibility({
      all: false,
      layerList: true,
      zoomControl: true,
      scaleControl: true,
      mapTypeControl: true,
      fullscreenControl: false
  });
  
  // Set the default map's cursor to a 'crosshair'.
  map.style().set('cursor', 'crosshair');
  
  // Set the center and zoom level of the new map.
  map.setCenter(4, 0, 3);

  // return the map object
  return map;
};


// function to create and set up the drawing tools object for the map
var createDrawingTools = function(map){
  var drawingTools = map.drawingTools();
  drawingTools.setShown(false);
  while (drawingTools.layers().length() > 0) {
    var layer = drawingTools.layers().get(0);
    drawingTools.layers().remove(layer);
  }
  return drawingTools;
};


// function for adding a reference to a panel
// Function to create reference panel.
function addRefToPanel(panel) {
    var referenceHeader = ui.Label({
        value: 'Data source:',
        style: {
            color: 'black',
            textAlign: 'center'
        },
    });
    var dataUseStatement = ui.Label({
        value: 'All data are free to use and redistribute, with proper ' +
               'attribution under the MIT license and with citation ' +
               'of the original scientific paper (URL; DOI).',
        style: {
            color: 'black',
            fontSize: '.7vw',
        },
    });
    var reference = ui.Label({
        value: 'PAPER TITLE HERE',
        style: {
            color: 'black',
            textAlign: 'center'
        },
        targetUrl: 'URL HERE'
    });

    // Add reference to the panel.
    panel.add(referenceHeader);
    panel.add(dataUseStatement);
    panel.add(reference);
}


// function to create layer selector
var createLayerSelector = function(map, EOFsForViz, asynch, mainPanel){

  // define info about each layer
  var dataInfo = {
    'EOFs': {
        name: 'LSP diversity',
        desc: 'false-color visualization of the major global diversity in LSP',
        img: EOFsForViz,
        vis: {
            min: 0,
            max: 1,
            opacity: 0.85
        },
        cbar: false,
    },
    'asynch': {
        name: 'LSP asynchrony',
        desc: 'spatial asynchrony in LSP (calculated within 100 km neighborhood radii surrounding each pixel)',
        img: asynch,
        vis: {
            min: 0,
            max: 0.0001465478,  // 99th percentile
            palette: asynchPalette,
            opacity: 0.85
        },
        cbar: true, 
    }
  };
  // create a layer selector, using that data info
  var items = [];
  Object.keys(dataInfo).forEach(function(key) {
    items.push({value: key, label: dataInfo[key].name});
  });
  items.push({value: 'none', label: 'Remove all'});

  var selector = ui.Select({
      items: items,
      value: items[0].value,
      style: {margin: '8px 0px'}
  });

  // create a redrawing function, to be called when the user changes the selected layer
  function redrawLayer(layer) {
    // Fetch the info that corresponds to the selected layer.
    var info = dataInfo[layer];

    // reset map layers
    map.layers().reset();
    map.remove(asynchLegend);
    
    // if layer is none, reset map layers
    if (layer == 'none') {
        map.layers().reset();
        map.remove(asynchLegend);
    } else {
    // add chosen layer to map
        var visImg = info.img.visualize(info.vis);
        map.remove(asynchLegend);
        map.addLayer(visImg, {}, layer);
        // add colorbar, if needed
        if (info.cbar){
          map.add(asynchLegend);
        }
    }
  }
  // register the redrawLayer function as a callback on the layer selector
  selector.onChange(redrawLayer);
  // add a header and then the layer selector to the main panel
  mainPanel.add(ui.Label('Layer selection:'));
  mainPanel.add(selector);
  // start off with the EOFs displayed
  redrawLayer('EOFs');
};


// function to create all the main panel's sub-panels
var createMainSubPanels = function(map, results, mainPanel, drawingTools){
  // define the summary-info panel
  var introPanel = ui.Panel([
      ui.Label({
            value: 'Global diversity and asynchrony of ecosystem seasonality',
            style: {
                fontSize: '1.1vw',
                fontWeight: 'bold'
                }
              }),
      ui.Label({
            value: 'This app provides access to the main mapping ' +
                   'products from Terasaki Hart et al. 2024. ' +
                   '(URL; DOI), which maps global estimates of ' +
                   'the diversity and asynchrony of annual vegetation seasonality ' +
                   '(also known as \'land surface phenology\' (LSP)) ' +
                   'and explores their ecological drivers and evolutionary implications. ',
            style: {fontSize: '.7vw',
                    fontWeight: 'normal',
                   },
              }),
      ui.Label('Instructions for use:'),
      ui.Label({value: '1. Use the \'Layer selection\' dropdown menu to add or remove map layers, ' +
                   '(and the layer list within the map panel to control their transparency, colors, etc, if desired).',
                style: {fontSize: '.5vw',
                        fontWeight: 'normal',
                },
              }),
      ui.Label({value: '2. In the \'Location selection\' section, click the \'🖍 Draw points\' button, then point and click anywhere on the map ' +
                   'to drop markers at any two locations.',
                style: {fontSize: '.5vw',
                        fontWeight: 'normal',
                },
              }),
      ui.Label({value: '3. Wait for the line plots to render, providing a visualizaiton of both the raw ' +
                   'and modeled annual LSP patterns at both points (based on MODIS near-infrared ' +
                   'reflectance of vegetation (NIRv)).',
                style: {fontSize: '.5vw',
                        fontWeight: 'normal',
                },
              }),
      ui.Label({value: '4. You can click and drag the markers to move them to new locations ' +
                   'and dynamically update the line plots.',
                style: {fontSize: '.5vw',
                        fontWeight: 'normal',
                },
              }),
      ui.Label({value: '5. Alternatively, you can click the \'🗑 ️Clear all points\' button ' +
                   'then repeat the above process, to select new locations.',
                style: {fontSize: '.5vw',
                        fontWeight: 'normal',
                },
              }),
        ]);
  // define on-click actions for both contol panel buttons
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
  // create control panel and its buttons
  var symbol = {
  pointDraw: '🖍️',
  pointErase: '🗑️',
  };
  var controlPanel = ui.Panel([
      ui.Label('Location selection:'),
      ui.Button({
        label: symbol.pointDraw + ' Draw points',
        onClick: drawPoint,
        style: {stretch: 'horizontal'}
        }),
      ui.Button({
        label: symbol.pointErase + ' Clear all points',
        onClick: clearGeometries,
        style: {stretch: 'horizontal'}
        })
      ]);
  // get the EOFs and asynch maps
  var EOFsForViz = results.select(['b10', 'b11', 'b12']);
  var asynch = results.select(['b13']);
  // put everything on the main panel
  mainPanel.add(introPanel);
  createLayerSelector(map, EOFsForViz, asynch, mainPanel);
  mainPanel.add(controlPanel);
  // add reference to the panel
  addRefToPanel(mainPanel);
};


// function to create the panel to hold the line plots
var createLinePlotPanel = function(drawingTools, LSPRawAndFittedTimeSeries){
  
  // create the line-plot panel and subpanels
  var lineplotPanel = ui.Panel({style: {position: 'bottom-left', width: '40%'}});
  var p1chartPanel = ui.Panel({
    style:
        {position: 'top-center', shown: false}
  });
  lineplotPanel.add(p1chartPanel);
  var p2chartPanel = ui.Panel({
    style:
        {position: 'bottom-center', shown: false}
  });
  lineplotPanel.add(p2chartPanel);  
  
  function formatPointCoordsString(pt, ptNumString){
    var lonString = ee.String(ee.Number(ee.List(pt.coordinates()).get(0)).format('%.2f'));
    var latString = ee.String(ee.Number(ee.List(pt.coordinates()).get(1)).format('%.2f'));
    var headerString = ee.String('POINT ').cat(ee.String(ptNumString)).cat(ee.String(': ('));
    var coordsString = headerString.cat(lonString).cat(ee.String(', ')).cat(latString).cat(ee.String(')'));
    return coordsString.getInfo();
  }
  
  function chartTimeSeries() {
    // only execute if at least two points have been drawn
    // (i.e., there are at least 4 total values in the coordinates array)
    var point_coords_array_len = drawingTools.layers().get(0).getEeObject().coordinates().flatten().length().getInfo(); 
    if (point_coords_array_len >= 4){
      // Make the chart panel visible the first time a geometry is drawn.
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
      
      // Set the drawing mode back to null; turns drawing off.
      drawingTools.setShape(null);
    
      // Chart time series for the selected area of interest.
      var p1chart = ui.Chart.image.series(LSPRawAndFittedTimeSeries.select(['NIRv', 'fitted']), p1)
        .setOptions({title: formatPointCoordsString(p1, '1'),
                     hAxis: {'title': ee.String(p1.coordinates())}});
      var p2chart = ui.Chart.image.series(LSPRawAndFittedTimeSeries.select(['NIRv', 'fitted']), p2)
        .setOptions({title: formatPointCoordsString(p2, '2'),
                     hAxis: {'title': p2.coordinates()}});
      
      // Replace the existing chart in the chart panel with the new chart.
      p1chartPanel.widgets().reset([p1chart]);
      p2chartPanel.widgets().reset([p2chart]);
    }
  }
  drawingTools.onDraw(ui.util.debounce(chartTimeSeries, 500));
  drawingTools.onEdit(ui.util.debounce(chartTimeSeries, 500));
  return lineplotPanel;
};
  
  
// main function for creating the app
var createApp = function(){
    
  // load results image
  var results = ee.Image('users/drewhart/Terasaki_Hart_2024_LSP_main_map_results');
  
  // load 20-year NIRv dataset
  var NIRv = getNIRvDataForViz();
  
  // calculate the harmonic regression fitted values
  var NIRvRawAndFitted = calcHarmonicRegressionForViz(NIRv);
  //print('NIRv time series (both raw and fitted values)', NIRvRawAndFitted);
  
  // configure to display our data over top of Google Earth satellite-imagery mosaic, with placename labels               
  Map.setOptions('HYBRID');
  Map.style().set('cursor', 'crosshair');
  Map.setCenter(4, 0, 3);
  
  // create the map
  var map = createMap();
  
  // set up the drawing tools
  var drawingTools = createDrawingTools(map);
  
  // create the main panel
  var mainPanel = ui.Panel({layout: ui.Panel.Layout.flow('vertical'),
                        style: {width: '15%', fontSize: '1vw'}
                      });
  
  // create the lineplot panel
  var lineplotPanel = createLinePlotPanel(drawingTools, NIRvRawAndFitted);
  
  // add the main panel to the root panel
  ui.root.insert(1, mainPanel);
  // add the lineplot panel to the map
  map.add(lineplotPanel);

  // call the subpanel creation function
  createMainSubPanels(map, results, mainPanel, drawingTools);
  
};  
  



/////////////////////////////////////////////////////////
// CREATE THE APP
/////////////////////////////////////////////////////////

// call the main app-creation function
createApp();


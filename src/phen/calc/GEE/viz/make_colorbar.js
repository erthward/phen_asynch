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

exports.makeLegend = function(title, low, mid, high, palette, pos) {
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
      position: pos
    });
  return panel;
};

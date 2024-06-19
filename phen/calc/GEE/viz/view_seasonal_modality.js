var modality = ee.Image('users/drewhart/seasonalModeSIF');
Map.addLayer(modality, {min:1, max:2, palette:['blue', 'white', 'orange']}, 'seasonal modality');
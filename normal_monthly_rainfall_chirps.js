// ==============================
// 1️⃣ Philippines boundary
// ==============================
var philippines = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017")
  .filter(ee.Filter.eq('country_na', 'Philippines'));

var studyArea = philippines.geometry();

Map.centerObject(studyArea, 6);
Map.addLayer(
  ee.Image().paint(philippines, 0, 1.5),
  {palette: ['ffffff']}, 'Philippines Boundary'
);

// ==============================
// 2️⃣ CHIRPS Daily — aggregate to May monthly totals
// ==============================
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .select('precipitation');

function monthlySum(year, month) {
  var start = ee.Date.fromYMD(year, month, 1);
  var end   = start.advance(1, 'month');
  return chirps
    .filterDate(start, end)
    .sum()
    .clip(studyArea)
    .set('year', year)
    .set('month', month)
    .set('system:time_start', start.millis());
}

// ==============================
// 3️⃣ Build Month collection 1991–2020
// ==============================
var years = ee.List.sequence(1991, 2020);

var MayCollection = ee.ImageCollection(
  years.map(function(y) {
    return monthlySum(ee.Number(y).toInt(), 5); // change month when needed 
  })
);

// ==============================
// 4️⃣ Mean Month rainfall (mm) across 30 years
// ==============================
var meanRainfall = MayCollection.mean().rename('Mean_Rainfall_mm');

// Additional statistics
var minRainfall  = MayCollection.min().rename('Min_Rainfall_mm');
var maxRainfall  = MayCollection.max().rename('Max_Rainfall_mm');
var stdDevRainfall = MayCollection
  .reduce(ee.Reducer.stdDev())
  .rename('StdDev_Rainfall_mm');

// ==============================
// 5️⃣ Visualize — mean rainfall
//    Philippines Month range typically 10–400 mm
// ==============================
var rainfallVis = {
  min: 0,
  max: 400,
  palette: [
    'ffffff',  // 0 mm     — dry
    'c6e9f7',  // ~40 mm
    '74c8e8',  // ~80 mm
    '2196f3',  // ~120 mm
    '1565c0',  // ~160 mm
    '0d47a1',  // ~200 mm
    '1b5e20',  // ~250 mm  — moderately wet
    '388e3c',  // ~300 mm
    'f9a825',  // ~340 mm
    'e65100',  // ~370 mm
    'b71c1c'   // ~400+ mm — very high
  ]
};

Map.addLayer(meanRainfall, rainfallVis,
  'Mean May Rainfall 1991–2020 (mm)');
Map.addLayer(minRainfall, rainfallVis,
  'Min May Rainfall 1991–2020 (mm)', false);
Map.addLayer(maxRainfall, rainfallVis,
  'Max May Rainfall 1991–2020 (mm)', false);
Map.addLayer(stdDevRainfall, {min: 0, max: 150, palette: ['fffde7','ff6f00']},
  'StdDev May Rainfall 1991–2020 (mm)', false);

// ==============================
// 6️⃣ Classify mean rainfall into categories
// ==============================
var classified = ee.Image(0)
  .where(meanRainfall.lt(25),                              1)  // Very dry
  .where(meanRainfall.gte(25).and(meanRainfall.lt(75)),    2)  // Dry
  .where(meanRainfall.gte(75).and(meanRainfall.lt(150)),   3)  // Below normal
  .where(meanRainfall.gte(150).and(meanRainfall.lt(250)),  4)  // Normal
  .where(meanRainfall.gte(250).and(meanRainfall.lt(350)),  5)  // Above normal
  .where(meanRainfall.gte(350),                            6)  // Very wet
  .clip(studyArea)
  .rename('RainfallClass');

Map.addLayer(classified, {
  min: 1, max: 6,
  palette: ['b71c1c','f46d43','fdae61','a6d96a','1a9850','0d47a1']
}, 'Rainfall Classes — Mean May 1991–2020', false);

// ==============================
// 7️⃣ Legend Panel
// ==============================
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '10px 14px',
    backgroundColor: 'rgba(255,255,255,0.93)'
  }
});

legend.add(ui.Label({
  value: '🌧️ Mean May Rainfall | 1991–2020 | Philippines',
  style: {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}
}));
legend.add(ui.Label({
  value: 'Source: CHIRPS Daily aggregated to monthly | 30-year mean',
  style: {fontSize: '10px', color: '#888', margin: '0 0 8px 0'}
}));

// Continuous color bar
var colorBar = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select(0)
    .multiply((400 - 0) / 100.0).add(0)
    .visualize({min: 0, max: 400, palette: rainfallVis.palette}),
  params: {bbox: '0,0,1,0.1', dimensions: '200x20'},
  style: {stretch: 'horizontal', margin: '0 0 4px 0', maxHeight: '20px'}
});
legend.add(colorBar);

var colorBarLabels = ui.Panel({
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {margin: '0 0 8px 0'}
});
colorBarLabels.add(ui.Label('0 mm',   {fontSize: '10px', color: '#555'}));
colorBarLabels.add(ui.Label('200 mm', {fontSize: '10px', color: '#555',
  textAlign: 'center', stretch: 'horizontal'}));
colorBarLabels.add(ui.Label('400+ mm', {fontSize: '10px', color: '#555'}));
legend.add(colorBarLabels);

// Class labels
var classItems = [
  {color: 'b71c1c', label: 'Very dry      (< 25 mm)'},
  {color: 'f46d43', label: 'Dry           (25–75 mm)'},
  {color: 'fdae61', label: 'Below normal  (75–150 mm)'},
  {color: 'a6d96a', label: 'Normal        (150–250 mm)'},
  {color: '1a9850', label: 'Above normal  (250–350 mm)'},
  {color: '0d47a1', label: 'Very wet      (≥ 350 mm)'}
];

legend.add(ui.Label({
  value: 'Rainfall classes (toggle layer):',
  style: {fontSize: '10px', color: '#666', margin: '4px 0 4px 0'}
}));

classItems.forEach(function(c) {
  var row = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
  row.add(ui.Label({
    style: {
      backgroundColor: '#' + c.color,
      padding: '8px', margin: '2px 6px 2px 0',
      width: '16px', height: '16px', border: '1px solid #aaa'
    }
  }));
  row.add(ui.Label({value: c.label, style: {fontSize: '11px', fontFamily: 'monospace'}}));
  legend.add(row);
});

legend.add(ui.Label({
  value: 'Toggle layers in Layers panel (top-right)',
  style: {fontSize: '9px', color: '#aaa', margin: '6px 0 0 0'}
}));

Map.add(legend);

// ==============================
// 8️⃣ Zonal Statistics
// ==============================
// National mean rainfall
print('🌧️ National mean May rainfall (mm):', meanRainfall.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: studyArea,
  scale: 1000,
  bestEffort: true
}));

// Area per rainfall class (km²)
var areaImage = ee.Image.pixelArea().divide(1e6)
  .addBands(classified);

var areaStats = areaImage.reduceRegion({
  reducer: ee.Reducer.sum().group({
    groupField: 1,
    groupName: 'class'
  }),
  geometry: studyArea,
  scale: 5566,
  maxPixels: 1e13,
  bestEffort: true
});

print('📊 Area per rainfall class (km²):', areaStats);

// ==============================
// 9️⃣ Time series chart — national mean May rainfall per year
// ==============================
var chart = ui.Chart.image.series({
  imageCollection: MayCollection,
  region: studyArea,
  reducer: ee.Reducer.mean(),
  scale: 5566,
  xProperty: 'system:time_start'
})
.setChartType('ColumnChart')
.setOptions({
  title: 'National Mean May Rainfall 1991–2020 — Philippines',
  hAxis: {title: 'Year', format: 'yyyy'},
  vAxis: {
    title: 'Rainfall (mm)',
    viewWindow: {min: 0}
  },
  series: {0: {color: '1565c0'}},
  backgroundColor: 'white',
  legend: {position: 'none'},
  bar: {groupWidth: '80%'}
});

print(chart);

// ==============================
// 🔟 Export rasters
// ==============================
Export.image.toDrive({
  image: meanRainfall,
  description: 'MeanRainfall_Philippines_May_1991_2020',
  folder: 'GEE_Exports',
  fileNamePrefix: 'MeanRainfall_Philippines_May_1991_2020',
  region: studyArea,
  scale: 1000,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: classified,
  description: 'RainfallClass_Philippines_May_1991_2020',
  folder: 'GEE_Exports',
  fileNamePrefix: 'RainfallClass_Philippines_May_1991_2020',
  region: studyArea,
  scale: 1000,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// ==============================
// Export — StdDev May Rainfall
// ==============================

Export.image.toDrive({
  image: stdDevRainfall,
  description: 'RainfallStdDev_Philippines_May_1991_2020',
  folder: 'GEE_Exports',
  fileNamePrefix: 'RainfallStdDev_Philippines_May_1991_2020',
  region: studyArea,
  scale: 1000,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// ==============================
// Export — Min May Rainfall
// ==============================
Export.image.toDrive({
  image: minRainfall,
  description: 'MinRainfall_Philippines_May_1991_2020',
  folder: 'GEE_Exports',
  fileNamePrefix: 'MinRainfall_Philippines_May_1991_2020',
  region: studyArea,
  scale: 5566,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// ==============================
// Export — Max May Rainfall
// ==============================
Export.image.toDrive({
  image: maxRainfall,
  description: 'MaxRainfall_Philippines_May_1991_2020',
  folder: 'GEE_Exports',
  fileNamePrefix: 'MaxRainfall_Philippines_May_1991_2020',
  region: studyArea,
  scale: 5566,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
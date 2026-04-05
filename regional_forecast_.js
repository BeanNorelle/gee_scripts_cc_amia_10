// ============================================================
// Region X — Projected July Rainfall Analysis V2.2  
// CHIRPS Daily | Climatology 1991–2020 | FAO GAUL Level-2
// ============================================================
// by: Bean Sob
// 
//
// PH-Rainfall-Projection Engine
// -----------------------------
// Description:
//     Generates spatially explicit rainfall projections for the Philippine 
//     territory by downscaling PAGASA provincial forecasts onto CHIRPS 
//     climatological rasters (1991-2020 baseline).

// Processing Logic:
//     1. Align CHIRPS raster data with Philippine administrative boundaries.
//     2. Map PAGASA's 'Percent of Normal' forecast to specific province geometries.
//     3. Execute pixel-wise linear scaling to produce the final projected depth.

// Inputs: 
//     - CHIRPS Monthly Average Raster (ARM)
//     - PAGASA Provincial Forecast Data (PRM)
// Output: 
//     - Projected Rainfall Raster (PRF)
//
//
// ============================================================
// USER INPUT VARIABLES
// ============================================================

// PAGASA forecast % of normal for each province in Region X (Northern Mindanao)
var bukidnonPercent = 88.4;
var camiguinPercent =  94.4;
var lanaoDelNortePercent = 86.6;
var misamisOccidentalPercent = 82.5;
var misamisOrientalPercent = 83.7;

var selectedMonth = 7;//change this value to select different month (1-12)

var monthsList = ['January','February','March','April','May','June','July','August','September','October','November','December'];

var monthString = monthsList[selectedMonth -1];

// ============================================================
// 1️⃣ PHILIPPINES BOUNDARY (full-detail LSIB 2017)
// ============================================================
var philippines = ee.FeatureCollection("USDOS/LSIB/2017")
  .filter(ee.Filter.eq('COUNTRY_NA', 'Philippines'));

var philippinesGeom = philippines.geometry();

// ============================================================
// 2️⃣ REGION X PROVINCES (FAO GAUL Level-2)
// ============================================================
var allProvinces = ee.FeatureCollection("FAO/GAUL/2015/level2")
  .filter(ee.Filter.eq('ADM0_NAME', 'Philippines'));
print('All PH region names:', 
  ee.FeatureCollection("FAO/GAUL/2015/level2")
    .filter(ee.Filter.eq('ADM0_NAME', 'Philippines'))
    .aggregate_array('ADM1_NAME')
    .distinct()
    .sort()
);


// print('All ADM1 regions:', allProvinces.aggregate_array('ADM1_NAME').distinct().sort());

var provinces = allProvinces
  .filter(ee.Filter.eq('ADM1_NAME', 'Region X (Northern Mindanao)'));

print('✅ Province count (should be 5):', provinces.size());
print('✅ Province names:', provinces.aggregate_array('ADM2_NAME').sort());

var regionGeom = provinces.geometry();

Map.centerObject(regionGeom, 8);

// ============================================================
// 3️⃣ CHIRPS DAILY → Monthly Sum Helper
// ============================================================
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .select('precipitation');

function monthlySum(year, month) {
  var start = ee.Date.fromYMD(year, month, 1);
  var end   = start.advance(1, 'month');
  return chirps
    .filterDate(start, end)
    .sum()
    .clip(philippinesGeom)
    .set('year',              year)
    .set('month',             month)
    .set('system:time_start', start.millis());
}

// ============================================================
// 4️⃣ Selected Month CLIMATOLOGY (1991–2020)
// ============================================================
var years = ee.List.sequence(1991, 2020);

var monthCollection = ee.ImageCollection(
  years.map(function(y) { return monthlySum(y, selectedMonth); }) //{selectedMonth}  got to line 37 if need to change month
);

var meanRainfall   = monthCollection.mean().rename('Mean_Rainfall_mm');
var stdDevRainfall = monthCollection
  .reduce(ee.Reducer.stdDev())
  .select('precipitation_stdDev')   // explicit band select before rename
  .rename('StdDev_Rainfall_mm');
var minRainfall    = monthCollection.min().rename('Min_Rainfall_mm');
var maxRainfall    = monthCollection.max().rename('Max_Rainfall_mm');

// Coefficient of Variation (%) — inter-annual reliability
var cvRainfall = stdDevRainfall
  .divide(meanRainfall)
  .multiply(100)
  .rename('CV_pct');

// ============================================================
// 5️⃣ FORECAST TABLE — % of Normal (PAGASA-style)
// ============================================================
var forecastTable = ee.FeatureCollection([
  ee.Feature(null, {ADM2_NAME: 'Bukidnon',          forecast_pct: bukidnonPercent}),
  ee.Feature(null, {ADM2_NAME: 'Camiguin',           forecast_pct: camiguinPercent}),
  ee.Feature(null, {ADM2_NAME: 'Lanao del Norte',    forecast_pct: lanaoDelNortePercent}),
  ee.Feature(null, {ADM2_NAME: 'Misamis Occidental', forecast_pct: misamisOccidentalPercent}),
  ee.Feature(null, {ADM2_NAME: 'Misamis Oriental',   forecast_pct: misamisOrientalPercent})
]);

// ============================================================
// 6️⃣ JOIN FORECAST → PROVINCES (safe, no-crash)
// ============================================================
var joined = provinces.map(function(prov) {
  var name  = prov.get('ADM2_NAME');
  var match = forecastTable
    .filter(ee.Filter.eq('ADM2_NAME', name))
    .first();

  var forecastPct = ee.Algorithms.If(match, ee.Number(match.get('forecast_pct')), 100);
  var droughtRisk = ee.Number(
    ee.Algorithms.If(
      ee.Number(forecastPct).lt(75),  'High',
      ee.Algorithms.If(
        ee.Number(forecastPct).lt(85), 'Moderate',
        ee.Algorithms.If(
          ee.Number(forecastPct).lt(95), 'Low', 'None'
        )
      )
    )
  );

  return prov
    .set('forecast_pct', forecastPct)
    .set('drought_risk',  droughtRisk);
});

print('📋 Joined province table:', joined.select(['ADM2_NAME','forecast_pct','drought_risk']));

// ============================================================
// 7️⃣ RASTERIZE FORECAST
// ============================================================
var forecastImage = joined.reduceToImage({
  properties: ['forecast_pct'],
  reducer: ee.Reducer.first()
}).rename('Forecast_pct');

// ============================================================
// 8️⃣ DERIVED RASTERS
// ============================================================
var projectedRainfall = meanRainfall
  .multiply(forecastImage.divide(100))
  .rename('Projected_Rainfall_mm')
  .clip(regionGeom);

var absoluteAnomaly = projectedRainfall
  .subtract(meanRainfall.clip(regionGeom))
  .rename('Absolute_Anomaly_mm');

var percentAnomaly = projectedRainfall
  .subtract(meanRainfall.clip(regionGeom))
  .divide(meanRainfall.clip(regionGeom))
  .multiply(100)
  .rename('Percent_Anomaly_pct');

var zScore = projectedRainfall
  .subtract(meanRainfall.clip(regionGeom))
  .divide(stdDevRainfall.clip(regionGeom))
  .rename('Z_Score');

// ============================================================
// 9️⃣ DROUGHT RISK CLASSIFICATION
//    Based on % of normal thresholds
// ============================================================
var droughtClass = ee.Image(0)
  .where(forecastImage.lt(75),                              3)  // High risk
  .where(forecastImage.gte(75).and(forecastImage.lt(85)),   2)  // Moderate risk
  .where(forecastImage.gte(85).and(forecastImage.lt(95)),   1)  // Low risk
  .where(forecastImage.gte(95),                             0)  // No risk
  .clip(regionGeom)
  .selfMask()
  .rename('Drought_Risk_Class');

// ============================================================
// 🔟 VISUALIZATION PALETTES  (matched to legend image)
// ============================================================

// 13 stops evenly spread across 0–400 mm
// Breakpoints ~0, 33, 67, 80, 107, 133, 160, 200, 233, 267, 300, 340, 370, 400
var rainfallVis = {
  min: 0, max: 400,
  palette: [
    'bf2a2a',   //   0 mm — Very Dry (dark red)
    'e04030',   //  33 mm — Very Dry (red)
    'e8762a',   //  67 mm — Dry (orange)
    'f5a040',   //  80 mm — Dry (light orange)
    'f5cc70',   // 107 mm — Below Normal (amber-yellow)
    'f5e8a0',   // 133 mm — Below Normal (pale yellow)
    'e8f0a8',   // 160 mm — Moderate (yellow-green)
    'b8dfa0',   // 200 mm — Above Normal (soft green)
    '7ec880',   // 233 mm — Above Normal (medium green)
    '50b890',   // 267 mm — Above Normal (teal-green)
    '48a8b0',   // 300 mm — Above Normal (teal)
    '5090c8',   // 340 mm — Very Wet (steel blue)
    '2868a8'    // 400 mm — Very Wet (medium blue)
  ]
};

var anomalyVis = {
  min: -100, max: 0,
  palette: ['bf2a2a','e8762a','f5e8a0','ffffff']  // red → orange → pale yellow → white
};

var zScoreVis = {
  min: -2, max: 0,
  palette: ['bf2a2a','e04030','f5cc70','ffffff']  // dark red → red → yellow → white
};

// Drought risk — aligns with the "% of normal" thresholds
var droughtVis = {
  min: 0, max: 3,
  palette: [
    'e8f0a8',   // 0 = No risk    → Moderate/light-green tone
    'f5e8a0',   // 1 = Low risk   → Below Normal / pale yellow
    'f5a040',   // 2 = Moderate   → Dry / light orange
    'e04030'    // 3 = High risk  → Very Dry / red
  ]
};

// ============================================================
// 1️⃣1️⃣ MAP LAYERS
// ============================================================
// Background: national mean rainfall (context)
Map.addLayer(
  meanRainfall.clip(philippinesGeom),
  rainfallVis,
  '🗺️ National Mean `' + monthString + '` Rainfall (1991–2020)', false
);

// Region X layers
Map.addLayer(projectedRainfall,  rainfallVis,   '🌧️ Projected ' + monthString + ' Rainfall (mm)');
Map.addLayer(absoluteAnomaly,    anomalyVis,    '📉 Absolute Anomaly (mm)',          false);
Map.addLayer(percentAnomaly,
  {min: -20, max: 0, palette: ['b71c1c','ff7043','fff9c4','ffffff']},
  '📊 Percent Anomaly (%)', false);
Map.addLayer(zScore,             zScoreVis,     '📐 Z-Score (climatological)',        false);
Map.addLayer(droughtClass,       droughtVis,    '⚠️ Drought Risk Classification',     false);
Map.addLayer(cvRainfall.clip(regionGeom),
  {min: 0, max: 100, palette: ['e0f7fa','0097a7','004d5e']},
  '📈 Coefficient of Variation (%)', false);

Map.addLayer(
  ee.Image().paint(provinces, 0, 1),
  {palette: ['222222']},
  '🗂️ Province Boundaries'
);

// ============================================================
// 1️⃣2️⃣ ZONAL STATISTICS PER PROVINCE
// ============================================================
var stats = projectedRainfall.reduceRegions({
  collection: provinces,
  reducer: ee.Reducer.mean()
    .combine(ee.Reducer.min(),    '', true)
    .combine(ee.Reducer.max(),    '', true)
    .combine(ee.Reducer.stdDev(), '', true),
  scale: 5566,
  crs: 'EPSG:4326'
});

var statsWithAll = stats.map(function(f) {
  var projMean = ee.Number(f.get('mean'));
  var climMean = meanRainfall.clip(regionGeom).reduceRegion({
    reducer:  ee.Reducer.mean(),
    geometry: f.geometry(),
    scale:    5566,
    bestEffort: true
  }).get('Mean_Rainfall_mm');

  var anomaly  = projMean.subtract(ee.Number(climMean));
  var pctAnom  = anomaly.divide(ee.Number(climMean)).multiply(100);

  return f
    .set('clim_mean_mm', climMean)
    .set('anomaly_mm',   anomaly)
    .set('anomaly_pct',  pctAnom)
    .set('forecast_pct', f.get('forecast_pct'));
});

print('📊 Province-level statistics:', statsWithAll.select([
  'ADM2_NAME','forecast_pct','clim_mean_mm','mean',
  'anomaly_mm','anomaly_pct','min','max','stdDev'
]));

// ============================================================
// 1️⃣3️⃣ CHARTS
// ============================================================
// A) Climatology vs Projected — grouped bar chart per province
var forecastImageForChart = joined.reduceToImage({
  properties: ['forecast_pct'],
  reducer: ee.Reducer.first()
});

var climStack = meanRainfall.clip(regionGeom)
  .rename('Climatology_mm');

var projStack = projectedRainfall.rename('Projected_mm');

var combined = climStack.addBands(projStack);

var barChart = ui.Chart.image.byRegion({
  image:      combined,
  regions:    provinces,
  reducer:    ee.Reducer.mean(),
  scale:      5566,
  xProperty:  'ADM2_NAME'
})
.setChartType('ColumnChart')
.setOptions({
  title:  monthString + ' Rainfall — Climatology vs Projected | Region X',
  hAxis:  {title: 'Province'},
  vAxis:  {title: 'Rainfall (mm)', viewWindow: {min: 0}},
  series: {
    0: {color: '1565c0', label: 'Climatology (1991–2020 mean)'},
    1: {color: 'e65100', label: 'Projected (% of normal applied)'}
  },
  isStacked:       false,
  bar:             {groupWidth: '70%'},
  backgroundColor: 'white',
  legend:          {position: 'top'},
  chartArea:       {width: '75%'}
});
print(barChart);

// B) Time series 1991–2020 with trendline
var tsChart = ui.Chart.image.series({
  imageCollection: monthCollection,
  region:          regionGeom,
  reducer:         ee.Reducer.mean(),
  scale:           5566,
  xProperty:       'system:time_start'
})
.setChartType('LineChart')
.setOptions({
  title:     'Regional Mean ' + monthString + ' Rainfall 1991–2020 — Region X',
  hAxis:     {title: 'Year', format: 'yyyy'},
  vAxis:     {title: 'Rainfall (mm)', viewWindow: {min: 0}},
  series:    {0: {color: '1565c0', lineWidth: 2, pointSize: 4}},
  trendlines:{0: {color: 'e65100', lineWidth: 1.5, opacity: 0.7,
                  visibleInLegend: true, labelInLegend: 'Trend'}},
  backgroundColor: 'white',
  curveType: 'function',
  legend:    {position: 'bottom'}
});
print(tsChart);

// C) % of Normal by province — horizontal bar
var pctChart = ui.Chart.feature.byFeature({
  features:  joined,
  xProperty: 'ADM2_NAME',
  yProperties: ['forecast_pct']
})
.setChartType('BarChart')
.setOptions({
  title:   '% of Normal — ' + monthString + ' Forecast | Region X',
  hAxis:   {title: '% of Normal', viewWindow: {min: 0, max: 110},
            baseline: 100, baselineColor: 'green'},
  vAxis:   {title: 'Province'},
  series:  {0: {color: '1565c0'}},
  backgroundColor: 'white',
  legend: {position: 'none'}
});
print(pctChart);

// ============================================================
// 1️⃣4️⃣ LEGEND & UI PANEL
// ============================================================
var panel = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '10px 14px',
    backgroundColor: 'rgba(255,255,255,0.95)',
    width: '280px'
  }
});

panel.add(ui.Label({
  value: '🌧️ Projected ' + monthString + ' Rainfall',
  style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 2px 0'}
}));
panel.add(ui.Label({
  value: 'Region X — Philippines | Based on PAGASA Forecast',
  style: {fontSize: '10px', color: '#777', margin: '0 0 6px 0'}
}));
panel.add(ui.Label({
  value: 'Source: CHIRPS Daily | Climatology: 1991–2020',
  style: {fontSize: '10px', color: '#aaa', margin: '0 0 10px 0'}
}));

// Continuous color bar
panel.add(ui.Label('Rainfall (mm)', {fontSize: '11px', fontWeight: 'bold', margin: '0 0 2px 0'}));
var colorBar = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select(0)
    .multiply(400 / 100.0)
    .visualize({min: 0, max: 400, palette: rainfallVis.palette}),
  params: {bbox: '0,0,1,0.1', dimensions: '240x18'},
  style: {stretch: 'horizontal', margin: '0 0 2px 0', maxHeight: '18px'}
});
panel.add(colorBar);

var cbLabels = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style: {margin: '0 0 10px 0'}});
cbLabels.add(ui.Label('0',    {fontSize: '10px', color: '#555'}));
cbLabels.add(ui.Label('200 mm', {fontSize: '10px', color: '#555', textAlign: 'center', stretch: 'horizontal'}));
cbLabels.add(ui.Label('400+', {fontSize: '10px', color: '#555'}));
panel.add(cbLabels);

// Drought risk legend
panel.add(ui.Label('⚠️ Drought Risk (toggle layer)', {fontSize: '11px', fontWeight: 'bold', margin: '0 0 4px 0'}));

var riskItems = [
  {color: 'f8d7da', label: 'High risk    — < 75% of normal'},
  {color: 'ffd0a0', label: 'Moderate     — 75–85% of normal'},
  {color: 'fff3cd', label: 'Low risk     — 85–95% of normal'},
  {color: 'd4edda', label: 'No risk      — ≥ 95% of normal'}
];

riskItems.forEach(function(item) {
  var row = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style: {margin: '1px 0'}});
  row.add(ui.Label({style: {
    backgroundColor: '#' + item.color,
    padding: '7px', margin: '2px 6px 2px 0',
    width: '14px', height: '14px', border: '1px solid #bbb'
  }}));
  row.add(ui.Label({value: item.label, style: {fontSize: '10px', fontFamily: 'monospace'}}));
  panel.add(row);
});

// Forecast % per province
panel.add(ui.Label('📋 Forecast (% of Normal)', {fontSize: '11px', fontWeight: 'bold', margin: '10px 0 4px 0'}));

var forecastData = [
  {name: 'Bukidnon',          pct: bukidnonPercent},
  {name: 'Camiguin',          pct: camiguinPercent},
  {name: 'Lanao del Norte',   pct: lanaoDelNortePercent},
  {name: 'Misamis Occidental',pct: misamisOccidentalPercent},
  {name: 'Misamis Oriental',  pct: misamisOrientalPercent}
];


forecastData.forEach(function(d) {
  var row = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style: {margin: '1px 0'}});
  row.add(ui.Label({value: d.name, style: {fontSize: '10px', stretch: 'horizontal'}}));
  row.add(ui.Label({value: d.pct + '%', style: {
    fontSize: '10px', fontWeight: 'bold',
    color: d.pct < 85 ? '#c62828' : d.pct < 90 ? '#e65100' : '#2e7d32'
  }}));
  panel.add(row);
});

// Rainfall  legend
panel.add(ui.Label('⚠️ Rainfall Amount', {fontSize: '11px', fontWeight: 'bold', margin: '10px 0 4px 0'}));

var rainfallItems = [
  {color: 'bf2a2a', label: '1 - 10 mm — Very Dry (dark red)'},
  {color: 'e04030', label: '10 - 40 mm — Very Dry (red)'},
  {color: 'e8762a', label: '40 - 60 mm — Dry (orange)'},
  {color: 'f5a040', label: '60 - 80 mm — Dry (light orange)'},
  {color: 'f5cc70', label: '80 - 100 mm — Below Normal (amber-yellow)'},
  {color: 'f5e8a0', label: '100 - 120 mm — Below Normal (pale yellow)'},
  {color: 'e8f0a8', label: '120 - 160 mm — Moderate (yellow-green)'},
  {color: 'b8dfa0', label: '160 - 200 mm — Above Normal (soft green)'},
  {color: '7ec880', label: '200 - 250 mm — Above Normal (medium green)'},
  {color: '50b890', label: '250 - 300 mm — Above Normal (teal-green)'},
  {color: '48a8b0', label: '300 - 340 mm — Above Normal (teal)'},
  {color: '5090c8', label: '340 - 370 mm — Very Wet (steel blue)'},
  {color: '2868a8', label: '370 - 400 mm — Very Wet (medium blue)'}
];

rainfallItems.forEach(function(item) {
  var row = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style: {margin: '1px 0'}});
  row.add(ui.Label({style: {
    backgroundColor: '#' + item.color,
    padding: '7px', margin: '2px 6px 2px 0',
    width: '14px', height: '14px', border: '1px solid #bbb'
  }}));
  row.add(ui.Label({value: item.label, style: {fontSize: '10px', fontFamily: 'monospace'}}));
  panel.add(row);
});

panel.add(ui.Label({
  value: 'Toggle layers in Layers panel (top-right)',
  style: {fontSize: '9px', color: '#bbb', margin: '8px 0 0 0'}
}));

Map.add(panel);

// ============================================================
// 1️⃣5️⃣ INSPECTOR: click map to get pixel values
// ============================================================
Map.style().set('cursor', 'crosshair');

var inspectorPanel = ui.Panel({
  style: {position: 'top-right', padding: '8px 12px', backgroundColor: 'rgba(255,255,255,0.95)'}
});
inspectorPanel.add(ui.Label('🔍 Click map to inspect', {fontSize: '11px', color: '#555'}));
Map.add(inspectorPanel);

Map.onClick(function(coords) {
  inspectorPanel.clear();
  inspectorPanel.add(ui.Label('📍 Inspecting...', {fontSize: '11px'}));

  var point = ee.Geometry.Point([coords.lon, coords.lat]);
  var stack = projectedRainfall
    .addBands(meanRainfall.clip(regionGeom))
    .addBands(absoluteAnomaly)
    .addBands(percentAnomaly)
    .addBands(forecastImage);

  var values = stack.reduceRegion({
    reducer:  ee.Reducer.first(),
    geometry: point,
    scale:    5566
  });

  values.evaluate(function(result) {
    inspectorPanel.clear();
    inspectorPanel.add(ui.Label('🔍 Pixel Inspector', {fontWeight: 'bold', fontSize: '12px'}));
    inspectorPanel.add(ui.Label(
      'Lon: ' + coords.lon.toFixed(4) + ' | Lat: ' + coords.lat.toFixed(4),
      {fontSize: '10px', color: '#777', margin: '0 0 6px 0'}
    ));

    var fields = [
      {key: 'Projected_Rainfall_mm', label: 'Projected (mm)'},
      {key: 'Mean_Rainfall_mm',      label: 'Climatology (mm)'},
      {key: 'Absolute_Anomaly_mm',   label: 'Anomaly (mm)'},
      {key: 'Percent_Anomaly_pct',   label: 'Anomaly (%)'},
      {key: 'Forecast_pct',          label: 'Forecast (% of normal)'}
    ];

    fields.forEach(function(f) {
      var val = result[f.key];
      var display = (val !== null && val !== undefined)
        ? (Math.round(val * 10) / 10).toString()
        : 'outside region';
      var row = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
      row.add(ui.Label(f.label + ':', {fontSize: '10px', stretch: 'horizontal', color: '#555'}));
      row.add(ui.Label(display, {fontSize: '10px', fontWeight: 'bold'}));
      inspectorPanel.add(row);
    });
  });
});

// ============================================================
// 1️⃣6️⃣ EXPORTS — all at CHIRPS native scale (5566 m)
// ============================================================
var exportBase = {
  region:    regionGeom,
  scale:     5566,
  crs:       'EPSG:4326',
  maxPixels: 1e13,
  folder:    'GEE_Exports'
};

// var exportLayers = [
//   {image: projectedRainfall,             name: 'Projected_July_Rainfall_RegionX'},
//   {image: meanRainfall.clip(regionGeom), name: 'Clim_Mean_July_Rainfall_RegionX'},
//   {image: absoluteAnomaly,               name: 'Absolute_Anomaly_mm_RegionX'},
//   {image: percentAnomaly,                name: 'Percent_Anomaly_pct_RegionX'},
//   {image: zScore,                        name: 'ZScore_RegionX'},
//   {image: droughtClass,                  name: 'Drought_Risk_Class_RegionX'},
//   {image: cvRainfall.clip(regionGeom),   name: 'CV_Rainfall_pct_RegionX'}
// ];

// exportLayers.forEach(function(layer) {
//   Export.image.toDrive(Object.assign({}, exportBase, {
//     image:          layer.image,
//     description:    layer.name,
//     fileNamePrefix: layer.name
//   }));
// });

Export.image.toDrive({
  image: projectedRainfall,
  description: 'Projected_' + monthString + '_Rainfall_RegionX',
  folder: 'GEE_Exports',
  fileNamePrefix: 'Projected_' + monthString + '_Rainfall_RegionX',
  region: regionGeom,
  scale: 5566,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: meanRainfall.clip(philippinesGeom),
  description: 'National_Mean_' + monthString + '_Rainfall_PH',
  folder: 'GEE_Exports',
  fileNamePrefix: 'National_Mean_' + monthString + '_Rainfall_PH',
  region: philippinesGeom,
  scale: 5566,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Export.image.toDrive({
//   image: meanRainfall,
//   description: 'Projected_July_Rainfall_RegionX',
//   folder: 'GEE_Exports',
//   fileNamePrefix: 'Projected_July_Rainfall_RegionX',
//   region: philippinesGeom,
//   scale: 5566,
//   crs: 'EPSG:4326',
//   maxPixels: 1e13
// });
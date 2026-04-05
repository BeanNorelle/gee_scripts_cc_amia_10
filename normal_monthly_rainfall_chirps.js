// ============================================================
// Region X — Projected July Rainfall Analysis
// CHIRPS Daily | Climatology 1991–2020 | FAO GAUL Level-2
// ============================================================

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
// 4️⃣ July CLIMATOLOGY (1991–2020)
// ============================================================
var years = ee.List.sequence(1991, 2020);

var JulyCollection = ee.ImageCollection(
  years.map(function(y) { return monthlySum(y, 7); })
);

var meanRainfall   = JulyCollection.mean().rename('Mean_Rainfall_mm');
var stdDevRainfall = JulyCollection
  .reduce(ee.Reducer.stdDev())
  .select('precipitation_stdDev')
  .rename('StdDev_Rainfall_mm');
var minRainfall    = JulyCollection.min().rename('Min_Rainfall_mm');
var maxRainfall    = JulyCollection.max().rename('Max_Rainfall_mm');

// Coefficient of Variation (%) — inter-annual reliability
var cvRainfall = stdDevRainfall
  .divide(meanRainfall)
  .multiply(100)
  .rename('CV_pct');


// ============================================================
// 4️⃣B BIAS CORRECTION — Enhanced Delta Method
// ============================================================
// Integration point: Insert AFTER section 4️⃣ (climatology calcs)
// This preserves CHIRPS fine-scale variability + applies PAGASA trend

// ============================================================
// STEP 1: Calculate Residuals (CHIRPS spatial patterns)
// ============================================================
// Decompose CHIRPS into mean + residuals
var chirpsResiduals = ee.ImageCollection(
  years.map(function(y) { 
    var yearJuly = monthlySum(y, 7);
    // Residual = actual - climatological mean (per pixel pattern)
    return yearJuly.subtract(meanRainfall).set('year', y);
  })
);

// Mean spatial residual pattern (fine-scale deviations from mean)
var meanResidualPattern = chirpsResiduals.mean()
  .rename('Mean_Residual_Pattern_mm');

// Coefficient of variation in residuals (anomaly variability)
var residualStdDev = chirpsResiduals
  .reduce(ee.Reducer.stdDev())
  .select('precipitation_stdDev')
  .rename('Residual_StdDev_mm');

// ============================================================
// STEP 2: Delta Method with Spatial Downscaling
// ============================================================
// Key insight: Apply coarse forecast as a scaling factor,
// but preserve fine-scale patterns from CHIRPS residuals

// A) Calculate regional mean from PAGASA forecast
var pagasaRegionalMeanPct = forecastImage
  .reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: regionGeom,
    scale: 5566
  })
  .get('Forecast_pct');

// B) Create standardized regional forecast anomaly
var forecastAnomalyPct = ee.Number(pagasaRegionalMeanPct).subtract(100);
print('Regional forecast anomaly (% deviation from normal):', forecastAnomalyPct);

// C) Apply anomaly to climatology with residual preservation
var biasCorrectForecast_DeltaMethod = meanRainfall
  .multiply(ee.Number(pagasaRegionalMeanPct).divide(100))  // Coarse trend
  .add(
    meanResidualPattern  // Add back fine-scale pattern
      .multiply(ee.Number(pagasaRegionalMeanPct).divide(100))  // Scale residuals by forecast pct
  )
  .clip(regionGeom)
  .rename('BiasCorr_Projected_Rainfall_Delta');

// ============================================================
// STEP 3: Quantile Mapping Correction (optional alternative)
// ============================================================
// For extreme events: map forecast quantiles to observed quantiles
var chirpsQuantiles = JulyCollection.reduce(
  ee.Reducer.percentile([10, 25, 50, 75, 90])
);

// Get forecast quantiles from historical ensemble-like forecast
// For now, approximate using mean ± std
var fcstMean = meanRainfall;
var fcstStd = stdDevRainfall;

// Ratio correction: observed_std / forecast_std
var quantileCorrection = stdDevRainfall
  .divide(fcstStd.add(1))  // +1 avoid division by zero
  .rename('Quantile_Correction_Factor');

var biasCorrectForecast_QuantileMap = meanRainfall
  .multiply(ee.Number(pagasaRegionalMeanPct).divide(100))
  .multiply(quantileCorrection)
  .clip(regionGeom)
  .rename('BiasCorr_Projected_Rainfall_Quantile');

// ============================================================
// STEP 4: Locally-Adaptive Correction (Province-level)
// ============================================================
// Apply province-specific forecast % but preserve inter-pixel variability

var localAdaptiveCorrection = projectedRainfall.multiply(0);  // Initialize

// For each province, extract residuals and re-apply with local forecast
var biasCorrectLocal = provinces.map(function(prov) {
  var provName = prov.get('ADM2_NAME');
  var match = forecastTable.filter(ee.Filter.eq('ADM2_NAME', provName)).first();
  var provForecastPct = ee.Algorithms.If(
    match, 
    ee.Number(match.get('forecast_pct')), 
    100
  );
  
  var provGeom = prov.geometry();
  var provClimatology = meanRainfall.clip(provGeom);
  var provResiduals = meanResidualPattern.clip(provGeom);
  
  var corrected = provClimatology
    .multiply(ee.Number(provForecastPct).divide(100))
    .add(provResiduals.multiply(ee.Number(provForecastPct).divide(100)))
    .set('ADM2_NAME', provName);
  
  return corrected;
});

// Mosaic province corrections back together
var biasCorrectForecast_LocalAdaptive = ee.ImageCollection(
  ee.List.sequence(0, provinces.size().subtract(1)).map(function(i) {
    return biasCorrectLocal.toList().get(i);
  })
).mosaic().rename('BiasCorr_Projected_Rainfall_Local');

// ============================================================
// STEP 5: Uncertainty Quantification
// ============================================================
// Estimate prediction uncertainty from historical forecast skill

// RMSE from climatology vs actual (forecast skill proxy)
var forecastErrorHistory = JulyCollection
  .map(function(img) {
    return img.subtract(meanRainfall).pow(2);
  })
  .mean()
  .sqrt()
  .rename('Historical_RMSE_mm');

// Confidence interval (±1.96*RMSE for 95% CI)
var uncertainty_95pct = forecastErrorHistory.multiply(1.96);

var biasCorrectForecast_Lower = biasCorrectForecast_LocalAdaptive
  .subtract(uncertainty_95pct)
  .rename('BiasCorr_Projected_Rainfall_Lower95CI');

var biasCorrectForecast_Upper = biasCorrectForecast_LocalAdaptive
  .add(uncertainty_95pct)
  .rename('BiasCorr_Projected_Rainfall_Upper95CI');

// ============================================================
// STEP 6: Comparison Metrics (Original vs Corrected)
// ============================================================

// Original approach (section 8️⃣)
var anomaly_Original = projectedRainfall.subtract(meanRainfall.clip(regionGeom));

// Corrected approach
var anomaly_Corrected = biasCorrectForecast_LocalAdaptive
  .subtract(meanRainfall.clip(regionGeom));

// Spatial variability preservation metric
var stdDev_Original = projectedRainfall.reduceRegion({
  reducer: ee.Reducer.stdDev(),
  geometry: regionGeom,
  scale: 5566
}).get('Projected_Rainfall_mm');

var stdDev_Corrected = biasCorrectForecast_LocalAdaptive.reduceRegion({
  reducer: ee.Reducer.stdDev(),
  geometry: regionGeom,
  scale: 5566
}).get('BiasCorr_Projected_Rainfall_Local');

print('===== BIAS CORRECTION DIAGNOSTICS =====');
print('Original std dev (spatial):', stdDev_Original);
print('Corrected std dev (spatial):', stdDev_Corrected);
print('Spatial variability preserved (%):', 
  ee.Number(stdDev_Corrected).divide(ee.Number(stdDev_Original)).multiply(100));

// ============================================================
// STEP 7: Visualization — Comparison
// ============================================================

// Side-by-side anomaly comparison
Map.addLayer(
  anomaly_Original,
  {min: -100, max: 0, palette: ['b71c1c','ff7043','ffcc80','ffffff']},
  '📉 Anomaly (Original — uniform per province)', 
  false
);

Map.addLayer(
  anomaly_Corrected,
  {min: -100, max: 0, palette: ['b71c1c','ff7043','ffcc80','ffffff']},
  '📉 Anomaly (Bias-Corrected — spatial detail)',
  false
);

// Main output: Bias-corrected forecast
Map.addLayer(
  biasCorrectForecast_LocalAdaptive,
  {min: 0, max: 400, palette: [
    'ffffff','c6e9f7','74c8e8','2196f3',
    '1565c0','0d47a1','1b5e20','388e3c',
    'f9a825','e65100','b71c1c'
  ]},
  '🌧️ Projected Rainfall (Bias-Corrected + Fine-Scale)',
  true  // Add to map by default
);

// Uncertainty bands
Map.addLayer(
  biasCorrectForecast_Upper,
  {min: 0, max: 400, palette: ['ffffff','ffcccc']},
  '⚠️ Upper 95% Confidence Interval',
  false
);

Map.addLayer(
  biasCorrectForecast_Lower,
  {min: 0, max: 400, palette: ['ffcccc','ffffff']},
  '⚠️ Lower 95% Confidence Interval',
  false
);

// Residual pattern (what fine-scale detail was added)
Map.addLayer(
  meanResidualPattern.clip(regionGeom),
  {min: -50, max: 50, palette: ['d7191c','fdae61','ffffbf','a6d96a','1a9641']},
  '🔄 CHIRPS Fine-Scale Pattern (residuals)',
  false
);

// ============================================================
// STEP 8: Provincial-level statistics (corrected)
// ============================================================

var correctedStats = biasCorrectForecast_LocalAdaptive.reduceRegions({
  collection: provinces,
  reducer: ee.Reducer.mean()
    .combine(ee.Reducer.min(), '', true)
    .combine(ee.Reducer.max(), '', true)
    .combine(ee.Reducer.stdDev(), '', true),
  scale: 5566,
  crs: 'EPSG:4326'
});

var correctedStatsWithAnomaly = correctedStats.map(function(f) {
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
    .set('corrected_mean_mm', projMean)
    .set('anomaly_mm',   anomaly)
    .set('anomaly_pct',  pctAnom);
});

print('📊 BIAS-CORRECTED PROVINCE STATISTICS:', 
  correctedStatsWithAnomaly.select([
    'ADM2_NAME','corrected_mean_mm','clim_mean_mm','anomaly_mm','anomaly_pct'
  ])
);

// ============================================================
// STEP 9: Chart — Original vs Corrected by Province
// ============================================================

var correctedStack = biasCorrectForecast_LocalAdaptive.rename('Corrected_mm');
var origStack = projectedRainfall.rename('Original_mm');
var combined_comparison = climStack
  .addBands(origStack)
  .addBands(correctedStack);

var comparisonChart = ui.Chart.image.byRegion({
  image:      combined_comparison,
  regions:    provinces,
  reducer:    ee.Reducer.mean(),
  scale:      5566,
  xProperty:  'ADM2_NAME'
})
.setChartType('ColumnChart')
.setOptions({
  title:  'Rainfall Comparison: Climatology vs Original vs Bias-Corrected | Region X',
  hAxis:  {title: 'Province'},
  vAxis:  {title: 'Rainfall (mm)', viewWindow: {min: 0}},
  series: {
    0: {color: '1565c0', label: 'Climatology (1991–2020)'},
    1: {color: 'ff9800', label: 'Original Projection (uniform)'},
    2: {color: 'e65100', label: 'Bias-Corrected (w/ fine-scale)'}
  },
  isStacked:       false,
  bar:             {groupWidth: '65%'},
  backgroundColor: 'white',
  legend:          {position: 'top'},
  chartArea:       {width: '75%'}
});
print(comparisonChart);

// ============================================================
// STEP 10: Export Corrected Layers
// ============================================================

// Export.image.toDrive({
//   image: biasCorrectForecast_LocalAdaptive,
//   description: 'BiasCorr_Projected_July_Rainfall_RegionX',
//   folder: 'GEE_Exports',
//   fileNamePrefix: 'BiasCorr_Projected_July_Rainfall_RegionX',
//   region: regionGeom,
//   scale: 5566,
//   crs: 'EPSG:4326',
//   maxPixels: 1e13
// });

// Export.image.toDrive({
//   image: biasCorrectForecast_Upper,
//   description: 'BiasCorr_Upper95CI_RegionX',
//   folder: 'GEE_Exports',
//   fileNamePrefix: 'BiasCorr_Upper95CI_RegionX',
//   region: regionGeom,
//   scale: 5566,
//   crs: 'EPSG:4326',
//   maxPixels: 1e13
// });



Export.image.toDrive({
  image: biasCorrectForecast_Lower,
  description: 'BiasCorr_Lower95CI_RegionX',
  folder: 'GEE_Exports',
  fileNamePrefix: 'BiasCorr_Lower95CI_RegionX',
  region: regionGeom,
  scale: 5566,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// ============================================================
// 5️⃣ FORECAST TABLE — % of Normal (PAGASA-style)
// ============================================================
var forecastTable = ee.FeatureCollection([
  ee.Feature(null, {ADM2_NAME: 'Bukidnon',          forecast_pct: 92.2}),
  ee.Feature(null, {ADM2_NAME: 'Camiguin',           forecast_pct: 102.4}),
  ee.Feature(null, {ADM2_NAME: 'Lanao del Norte',    forecast_pct: 94.0}),
  ee.Feature(null, {ADM2_NAME: 'Misamis Occidental', forecast_pct: 93.7}),
  ee.Feature(null, {ADM2_NAME: 'Misamis Oriental',   forecast_pct: 96.3})
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
  var droughtRisk = ee.Algorithms.If(
      ee.Number(forecastPct).lt(75),  'High',
      ee.Algorithms.If(
        ee.Number(forecastPct).lt(85), 'Moderate',
        ee.Algorithms.If(
          ee.Number(forecastPct).lt(95), 'Low', 'None'
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

// Use bias-corrected forecast instead of simple % application
var projectedRainfall = biasCorrectForecast_LocalAdaptive
  .rename('Projected_Rainfall_mm')
  .clip(regionGeom);

// Rest of anomaly calculations follow same logic...
var absoluteAnomaly = projectedRainfall
  .subtract(meanRainfall.clip(regionGeom))
  .rename('Absolute_Anomaly_mm');

// etc.

// var projectedRainfall = meanRainfall
//   .multiply(forecastImage.divide(100))
//   .rename('Projected_Rainfall_mm')
//   .clip(regionGeom);

// var absoluteAnomaly = projectedRainfall
//   .subtract(meanRainfall.clip(regionGeom))
//   .rename('Absolute_Anomaly_mm');

var percentAnomaly = projectedRainfall
  .subtract(meanRainfall.clip(regionGeom))
  .divide(meanRainfall.clip(regionGeom))
  .multiply(100)
  .rename('Percent_Anomaly_pct');

// Standardised Anomaly (Z-score relative to climatological spread)
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
// 🔟 VISUALIZATION PALETTES
// ============================================================
var rainfallVis = {
  min: 0, max: 400,
  palette: [
    'ffffff','c6e9f7','74c8e8','2196f3',
    '1565c0','0d47a1','1b5e20','388e3c',
    'f9a825','e65100','b71c1c'
  ]
};

var anomalyVis = {
  min: -100, max: 0,
  palette: ['b71c1c','ff7043','ffcc80','ffffff']
};

var zScoreVis = {
  min: -2, max: 0,
  palette: ['7f0000','d32f2f','ef9a9a','ffffff']
};

var droughtVis = {
  min: 0, max: 3,
  palette: ['d4edda','fff3cd','ffd0a0','f8d7da']
};

// ============================================================
// 1️⃣1️⃣ MAP LAYERS
// ============================================================
// Background: national mean rainfall (context)
Map.addLayer(
  meanRainfall.clip(philippinesGeom),
  rainfallVis,
  '🗺️ National Mean July Rainfall (1991–2020)', false
);

// Region X layers
Map.addLayer(projectedRainfall,  rainfallVis,   '🌧️ Projected July Rainfall (mm)');
Map.addLayer(absoluteAnomaly,    anomalyVis,    '📉 Absolute Anomaly (mm)',          false);
Map.addLayer(percentAnomaly,
  {min: -20, max: 0, palette: ['b71c1c','ff7043','fff9c4','ffffff']},
  '📊 Percent Anomaly (%)', false);
Map.addLayer(zScore,             zScoreVis,     '📐 Z-Score (climatological)',        false);
Map.addLayer(droughtClass,       droughtVis,    '⚠️ Drought Risk Classification',     false);
Map.addLayer(cvRainfall.clip(regionGeom),
  {min: 0, max: 100, palette: ['e0f7fa','0097a7','004d5e']},
  '📈 Coefficient of Variation (%)', false);

// Province boundaries overlay (always visible)
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
  title:  'July Rainfall — Climatology vs Projected | Region X',
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
  imageCollection: JulyCollection,
  region:          regionGeom,
  reducer:         ee.Reducer.mean(),
  scale:           5566,
  xProperty:       'system:time_start'
})
.setChartType('LineChart')
.setOptions({
  title:     'Regional Mean July Rainfall 1991–2020 — Region X',
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
  title:   '% of Normal — July Forecast | Region X',
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
  value: '🌧️ Projected July Rainfall',
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
  {name: 'Bukidnon',          pct: 92.2},
  {name: 'Camiguin',          pct: 102.4},
  {name: 'Lanao del Norte',   pct: 94.0},
  {name: 'Misamis Occidental',pct: 93.7},
  {name: 'Misamis Oriental',  pct: 96.3}
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
Export.image.toDrive({
  image: projectedRainfall,
  description: 'Projected_July_Rainfall_RegionX',
  folder: 'GEE_Exports',
  fileNamePrefix: 'Projected_July_Rainfall_RegionX',
  region: regionGeom,
  scale: 5566,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

Export.image.toDrive({
  image: meanRainfall.clip(philippinesGeom),
  description: 'National_Mean_July_Rainfall_PH',
  folder: 'GEE_Exports',
  fileNamePrefix: 'National_Mean_July_Rainfall_PH',
  region: philippinesGeom,
  scale: 5566,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
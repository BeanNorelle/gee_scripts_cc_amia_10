// ============================================================
// QUANTILE MAPPING BIAS CORRECTION
// Observed  : CHIRPS Daily  — April climatology 1991–2020
// Hindcast  : CFSv2 FOR6H  — April 2011–2020
// Target    : Region X (Northern Mindanao), Philippines
// ============================================================

// ============================================================
// CONFIG — all tuneable parameters in one place
// ============================================================
var CFG = {
  obsYearStart:  1991,
  obsYearEnd:    2020,
  hindYearStart: 2011,   // CFSv2 available from 2011
  hindYearEnd:   2023,
  targetMonth:   4,      // April
  nQuantiles:    70,     // increase to 50–100 for smoother mapping
  scale:         5566,   // CHIRPS native resolution (m)
  crs:           'EPSG:4326',
  exportFolder:  'GEE_Exports',
  region:        'Region X (Northern Mindanao)'
};

// Forecast % of normal per province (PAGASA-issued)
var FORECAST_TABLE = [
  {name: 'Bukidnon',          pct: 83.1},
  {name: 'Camiguin',          pct: 87.6},
  {name: 'Lanao del Norte',   pct: 88.5},
  {name: 'Misamis Occidental',pct: 88.0},
  {name: 'Misamis Oriental',  pct: 85.1}
];

// ============================================================
// 1️⃣ BOUNDARIES
// ============================================================
var provinces = ee.FeatureCollection('FAO/GAUL/2015/level2')
  .filter(ee.Filter.eq('ADM0_NAME', 'Philippines'))
  .filter(ee.Filter.eq('ADM1_NAME', CFG.region));

print('✅ Province count (expect 5):', provinces.size());
print('✅ Province names:', provinces.aggregate_array('ADM2_NAME').sort());

var regionGeom = provinces.geometry();
Map.centerObject(regionGeom, 8);

// ============================================================
// 2️⃣ CHIRPS OBSERVED COLLECTION
//    Band renamed to 'rain' early → consistent naming throughout
// ============================================================
var chirpsDailyRaw = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
  .select(['precipitation'], ['rain']);  // rename on load

function chirpsMonthlySum(year, month) {
  var start = ee.Date.fromYMD(year, month, 1);
  return chirpsDailyRaw
    .filterDate(start, start.advance(1, 'month'))
    .sum()
    .clip(regionGeom)
    .set('year', year)
    .set('system:time_start', start.millis());
}

var obsCollection = ee.ImageCollection(
  ee.List.sequence(CFG.obsYearStart, CFG.obsYearEnd)
    .map(function(y) { return chirpsMonthlySum(y, CFG.targetMonth); })
);

var climMean = obsCollection.mean().rename('clim_mean_mm');

print('📦 CHIRPS obs count (expect 30):', obsCollection.size());

// ============================================================
// 3️⃣ CFSv2 HINDCAST COLLECTION
//    FIX: filter directly by month — no day-loop needed
//    FIX: band renamed to 'rain' immediately → matches CHIRPS naming
//    Rate (kg m⁻² s⁻¹) × 21600 s per 6-hr step = mm per image
// ============================================================
var cfsRaw = ee.ImageCollection('NOAA/CFSV2/FOR6H')
  .select(
    ['Precipitation_rate_surface_6_Hour_Average'],
    ['rain']   // rename on load — fixes band-name mismatch in QM function
  )
  .map(function(img) {
    return img.multiply(21600)   // rate → mm per 6-hr window
              .copyProperties(img, ['system:time_start']);
  });

function cfsMonthlySum(year, month) {
  var start = ee.Date.fromYMD(year, month, 1);
  return cfsRaw
    .filterDate(start, start.advance(1, 'month'))
    .sum()                           // sum all 6-hr images in April
    .clip(regionGeom)
    .set('year', year)
    .set('system:time_start', start.millis());
}

var hindcastCollection = ee.ImageCollection(
  ee.List.sequence(CFG.hindYearStart, CFG.hindYearEnd)
    .map(function(y) { return cfsMonthlySum(y, CFG.targetMonth); })
);

print('📦 CFSv2 hindcast count (expect 10):', hindcastCollection.size());

// ============================================================
// 4️⃣ RAW FORECAST IMAGE
//    % of normal applied per province to the 30-yr climatology
// ============================================================
var forecastFC = ee.FeatureCollection(
  FORECAST_TABLE.map(function(d) {
    return ee.Feature(null, {ADM2_NAME: d.name, forecast_pct: d.pct});
  })
);

var forecastPctImage = provinces.map(function(prov) {
  var match = forecastFC
    .filter(ee.Filter.eq('ADM2_NAME', prov.get('ADM2_NAME')))
    .first();
  var pct = ee.Algorithms.If(match, match.get('forecast_pct'), 100);
  return prov.set('forecast_pct', pct);
}).reduceToImage({
  properties: ['forecast_pct'],
  reducer:    ee.Reducer.first()
});

var rawForecast = climMean
  .multiply(forecastPctImage.divide(100))
  .rename('rain');   // keep band name consistent

// ============================================================
// 5️⃣ PERCENTILE BREAKPOINTS
//    Both collections share 'rain' band → identical band names
//    in the outputs: rain_p0, rain_p5, … rain_p100
// ============================================================
var percentileSteps = ee.List.sequence(0, 100, 100 / CFG.nQuantiles);

function buildPercentiles(collection) {
  return collection.reduce(
    ee.Reducer.percentile(
      percentileSteps.map(function(p) { return ee.Number(p).int(); })
    )
  ).clip(regionGeom);
}

var obsPercentiles  = buildPercentiles(obsCollection);
var fcstPercentiles = buildPercentiles(hindcastCollection);

print('🔢 Obs percentile bands:',  obsPercentiles.bandNames());
print('🔢 Fcst percentile bands:', fcstPercentiles.bandNames());

// ============================================================
// 6️⃣ QUANTILE MAPPING — piecewise linear interpolation
//    Both obs and forecast use 'rain_p*' naming → no mismatch
// ============================================================
function quantileMap(rawImage, fcstPerc, obsPerc, nSteps) {

  var corrected = obsPerc.select('rain_p0');   // floor = lowest obs quantile

  for (var i = 0; i < nSteps; i++) {
    var pLow  = Math.round(i       * (100 / nSteps));
    var pHigh = Math.round((i + 1) * (100 / nSteps));

    var fLo = fcstPerc.select('rain_p' + pLow);
    var fHi = fcstPerc.select('rain_p' + pHigh);
    var oLo = obsPerc.select('rain_p'  + pLow);
    var oHi = obsPerc.select('rain_p'  + pHigh);

    var inBracket   = rawImage.gte(fLo).and(rawImage.lt(fHi));
    var fRange      = fHi.subtract(fLo).max(0.001);  // guard div-by-zero
    var fraction    = rawImage.subtract(fLo).divide(fRange).clamp(0, 1);
    var mappedValue = oLo.add(fraction.multiply(oHi.subtract(oLo)));

    corrected = corrected.where(inBracket, mappedValue);
  }

  // Extrapolate above p100 using the observed ceiling
  corrected = corrected.where(
    rawImage.gte(fcstPerc.select('rain_p100')),
    obsPerc.select('rain_p100')
  );

  return corrected.rename('qm_corrected_mm');
}

var qmCorrected = quantileMap(
  rawForecast,
  fcstPercentiles,
  obsPercentiles,
  CFG.nQuantiles
).clip(regionGeom);

// ============================================================
// 7️⃣ VERIFICATION STATISTICS
// ============================================================
function regionalMean(image, label) {
  return image.reduceRegion({
    reducer:    ee.Reducer.mean(),
    geometry:   regionGeom,
    scale:      CFG.scale,
    bestEffort: true
  }).set('label', label);
}

print('📊 Regional mean — Climatology (mm):',  regionalMean(climMean,    'climatology'));
print('📊 Regional mean — Raw forecast (mm):',  regionalMean(rawForecast, 'raw'));
print('📊 Regional mean — QM corrected (mm):',  regionalMean(qmCorrected, 'qm_corrected'));

// Province-level bias table
var biasTable = qmCorrected.rename('qm')
  .addBands(rawForecast.rename('raw'))
  .addBands(climMean.rename('clim'))
  .reduceRegions({
    collection: provinces,
    reducer:    ee.Reducer.mean(),
    scale:      CFG.scale
  });

print('📋 Province bias table:', biasTable.select(['ADM2_NAME','qm','raw','clim']));

// ============================================================
// 8️⃣ CROSS-VALIDATION — apply QM to each hindcast year
//    Lets you see how well QM corrects known years
// ============================================================
var crossVal = ee.ImageCollection(
  ee.List.sequence(CFG.hindYearStart, CFG.hindYearEnd).map(function(y) {
    var rawYear = cfsMonthlySum(y, CFG.targetMonth);
    var qmYear  = quantileMap(rawYear, fcstPercentiles, obsPercentiles, CFG.nQuantiles)
                    .rename('qm');
    var obsYear = chirpsMonthlySum(y, CFG.targetMonth).rename('obs');

    return rawYear.rename('raw')
      .addBands(qmYear)
      .addBands(obsYear)
      .set('system:time_start', ee.Date.fromYMD(y, CFG.targetMonth, 1).millis());
  })
);

var crossValChart = ui.Chart.image.series({
  imageCollection: crossVal,
  region:          regionGeom,
  reducer:         ee.Reducer.mean(),
  scale:           CFG.scale,
  xProperty:       'system:time_start'
})
.setChartType('LineChart')
.setOptions({
  title:      'Cross-validation: CFSv2 raw vs QM-corrected vs CHIRPS observed',
  hAxis:      {title: 'Year', format: 'yyyy'},
  vAxis:      {title: 'April rainfall (mm)', viewWindow: {min: 0}},
  series: {
    0: {color: 'e65100', lineWidth: 1.5, lineDashStyle: [4, 2],
        pointSize: 3, visibleInLegend: true},
    1: {color: '1565c0', lineWidth: 2,
        pointSize: 4, visibleInLegend: true},
    2: {color: '2e7d32', lineWidth: 1.5,
        pointSize: 3, visibleInLegend: true}
  },
  trendlines: {
    0: {color: 'e65100', opacity: 0.25, lineWidth: 1},
    1: {color: '1565c0', opacity: 0.25, lineWidth: 1}
  },
  series: {
    0: {color: 'e65100', label: 'CFSv2 raw'},
    1: {color: '1565c0', label: 'QM corrected'},
    2: {color: '2e7d32', label: 'CHIRPS observed'}
  },
  backgroundColor: 'white',
  legend: {position: 'top'}
});

print(crossValChart);

// ============================================================
// 9️⃣ MAP LAYERS
// ============================================================
var rainfallVis = {
  min: 0, max: 400,
  palette: [
    'ffffff','c6e9f7','74c8e8','2196f3',
    '1565c0','0d47a1','1b5e20','388e3c',
    'f9a825','e65100','b71c1c'
  ]
};

var deltaVis = {
  min: -80, max: 80,
  palette: ['b71c1c','ff7043','fff9c4','aed6f1','1565c0']
};

Map.addLayer(climMean,    rainfallVis, '📘 Climatology mean (mm)');
Map.addLayer(rawForecast, rainfallVis, '🟠 Raw forecast (mm)');
Map.addLayer(qmCorrected, rainfallVis, '🟢 QM-corrected forecast (mm)');
Map.addLayer(
  qmCorrected.subtract(rawForecast.rename('qm_corrected_mm')).rename('delta'),
  deltaVis,
  '🔄 QM correction delta (mm)', false
);
Map.addLayer(
  hindcastCollection.mean().rename('cfs_mean'),
  rainfallVis,
  '🔵 CFSv2 hindcast mean (mm)', false
);
Map.addLayer(
  ee.Image().paint(provinces, 0, 1),
  {palette: ['222222']},
  '🗂️ Province boundaries'
);

// ============================================================
// 🔟 LEGEND
// ============================================================
var legend = ui.Panel({
  style: {
    position: 'bottom-left', padding: '10px 14px',
    backgroundColor: 'rgba(255,255,255,0.95)', width: '270px'
  }
});

legend.add(ui.Label('QM Bias Correction — April | Region X', {
  fontWeight: 'bold', fontSize: '13px', margin: '0 0 2px 0'
}));
legend.add(ui.Label(
  'Obs: CHIRPS ' + CFG.obsYearStart + '–' + CFG.obsYearEnd +
  ' | Hindcast: CFSv2 ' + CFG.hindYearStart + '–' + CFG.hindYearEnd,
  {fontSize: '10px', color: '#888', margin: '0 0 8px 0'}
));

[
  {color: '2e7d32', label: 'Climatology mean (CHIRPS)'},
  {color: 'e65100', label: 'Raw forecast (% of normal)'},
  {color: '1565c0', label: 'QM-corrected forecast'},
].forEach(function(item) {
  var row = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style: {margin: '2px 0'}});
  row.add(ui.Label({style: {
    backgroundColor: '#' + item.color, padding: '7px',
    margin: '2px 6px 2px 0', width: '14px', height: '14px',
    border: '1px solid #bbb'
  }}));
  row.add(ui.Label(item.label, {fontSize: '11px'}));
  legend.add(row);
});

legend.add(ui.Label(
  'Quantile steps: ' + CFG.nQuantiles + ' | Scale: ' + CFG.scale + ' m',
  {fontSize: '9px', color: '#bbb', margin: '8px 0 0 0'}
));

Map.add(legend);

// ============================================================
// 1️⃣1️⃣ EXPORTS
// ============================================================
var exportBase = {
  region: regionGeom, scale: CFG.scale,
  crs: CFG.crs, maxPixels: 1e13, folder: CFG.exportFolder
};

[
  {image: qmCorrected,          name: 'QM_Corrected_April_RegionX'},
  {image: rawForecast.rename('raw_forecast_mm'),
                                name: 'Raw_Forecast_April_RegionX'},
  {image: climMean,             name: 'Climatology_Mean_April_RegionX'},
  {image: obsPercentiles,       name: 'CHIRPS_Percentiles_April_RegionX'},
  {image: fcstPercentiles,      name: 'CFSv2_Percentiles_April_RegionX'},
].forEach(function(layer) {
  Export.image.toDrive(Object.assign({}, exportBase, {
    image:          layer.image,
    description:    layer.name,
    fileNamePrefix: layer.name
  }));
});
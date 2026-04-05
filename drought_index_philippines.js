// ============================================================
// DROUGHT RISK INDEX — Philippines | 1991–2020
// Indicators: SPI + NDVI Anomaly + Land Cover/Crop Exposure
// Risk = weighted composite of normalized indicators
// ============================================================

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

// ============================================================
// INDICATOR 1 — SPI (Standardized Precipitation Index)
// Represents rainfall deficit / drought frequency
// ============================================================

// ==============================
// 2️⃣ CHIRPS Daily → annual April monthly sums
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
    .set('system:time_start', start.millis());
}

// Build dry season months collection (Feb–May = peak dry season PH)
// We average SPI across dry season months for a robust drought signal
var dryMonths = [2, 3, 4, 5];
var years = ee.List.sequence(1991, 2020);

var drySeasonCollection = ee.ImageCollection(
  years.map(function(y) {
    var monthImages = dryMonths.map(function(m) {
      return monthlySum(ee.Number(y).toInt(), m);
    });
    // Mean precip across Feb-May for this year
    return ee.ImageCollection(monthImages).mean()
      .set('year', y)
      .set('system:time_start',
        ee.Date.fromYMD(ee.Number(y).toInt(), 3, 1).millis());
  })
);

// Baseline stats
var precipMean   = drySeasonCollection.mean();
var precipStdDev = drySeasonCollection.reduce(ee.Reducer.stdDev());

// Mean SPI across all years (pixel-wise mean of annual SPIs)
var spiCollection = ee.ImageCollection(
  years.map(function(y) {
    var precip = drySeasonCollection
      .filter(ee.Filter.eq('year', y))
      .first();
    return precip.subtract(precipMean)
      .divide(precipStdDev)
      .rename('SPI')
      .set('year', y);
  })
);

var meanSPI = spiCollection.mean().rename('Mean_SPI');

// Drought frequency — fraction of years with SPI < -1.0 (moderately dry)
var droughtFreq = spiCollection
  .map(function(img) { return img.lt(-1.0).rename('drought_yr'); })
  .sum()
  .divide(30)
  .rename('Drought_Frequency')
  .clip(studyArea);

print('✅ SPI computed');

// ============================================================
// INDICATOR 2 — NDVI Anomaly (Vegetation Stress)
// Negative anomaly = vegetation under stress = drought signal
// ============================================================

// ==============================
// 3️⃣ MODIS MOD13A3 — monthly NDVI (1km)
//    Dry season (Feb–May) mean NDVI anomaly
// ==============================

var modisNDVI = ee.ImageCollection('MODIS/061/MOD13A3')
  .select('NDVI')
  .filter(ee.Filter.calendarRange(2, 5, 'month'))
  .filterDate('2001-01-01', '2021-01-01')
  .filterBounds(studyArea);  // spatial filter without reprojection

// Scale NDVI — do NOT clip inside map() to avoid CRS conflict
var ndviScaled = modisNDVI.map(function(img) {
  return img.multiply(0.0001)
    .rename('NDVI')
    .copyProperties(img, ['system:time_start']);
});

// Compute baseline stats at native MODIS projection
var ndviMean   = ndviScaled.mean();
var ndviStdDev = ndviScaled.reduce(ee.Reducer.stdDev());

// Standardized anomaly per image, then mean across all years
var ndviAnomaly = ndviScaled
  .map(function(img) {
    return img.subtract(ndviMean)
      .divide(ndviStdDev)
      .rename('NDVI_Anomaly')
      .copyProperties(img, ['system:time_start']);
  })
  .mean()
  .rename('Mean_NDVI_Anomaly');

// Flip sign: negative anomaly = stress. Clip ONLY here at the end
var ndviStress = ndviAnomaly
  .multiply(-1)
  .rename('NDVI_Stress')
  .clip(studyArea);  // ✅ clip after all reductions, not inside map()

print('✅ NDVI anomaly computed');

// Normalize NDVI stress (-3 to 3 → 0 to 1)
var ndviNorm = ndviStress
  .subtract(-3)
  .divide(6)
  .clamp(0, 1)
  .rename('NDVI_norm');
// ============================================================
// INDICATOR 3 — Land Cover / Crop Exposure
// Cropland pixels = higher drought risk to food security
// ============================================================

// ==============================
// 4️⃣ ESA WorldCover v200 (2021) — cropland mask
// ==============================
var worldCover = ee.ImageCollection('ESA/WorldCover/v200')
  .first()
  .clip(studyArea);

// Class 40 = Cropland
// Class 10 = Tree cover (coconut/orchard proxy — also drought sensitive)
// Class 30 = Grassland (pasture — moderate sensitivity)
var cropExposure = ee.Image(0)
  .where(worldCover.eq(10), 0.5)   // Tree cover — moderate exposure
  .where(worldCover.eq(30), 0.6)   // Grassland  — moderate-high
  .where(worldCover.eq(40), 1.0)   // Cropland   — highest exposure
  .rename('Crop_Exposure')
  .clip(studyArea);

print('✅ Land cover exposure computed');

// ============================================================
// 5️⃣ NORMALIZE all indicators to 0–1 scale
//    0 = lowest risk, 1 = highest risk
// ============================================================
function normalize(image, bandName, minVal, maxVal) {
  return image.subtract(minVal)
    .divide(maxVal - minVal)
    .clamp(0, 1)
    .rename(bandName);
}

// SPI indicator: drought frequency (already 0–1, higher = more drought)
var spiNorm = droughtFreq.clamp(0, 1).rename('SPI_norm');

// NDVI stress: standardized anomaly flipped, clamp -3 to 3 → 0 to 1
var ndviNorm = normalize(ndviStress, 'NDVI_norm', -3, 3);

// Crop exposure: already 0–1
var cropNorm = cropExposure.rename('Crop_norm');

// ============================================================
// 6️⃣ WEIGHTED COMPOSITE DROUGHT RISK INDEX
//    Weights reflect contribution to agricultural drought risk
//    SPI:  40% — primary driver (rainfall deficit)
//    NDVI: 35% — vegetation response to drought
//    Crop: 25% — exposure/sensitivity of land use
// ============================================================
var droughtRisk = spiNorm.multiply(0.40)
  .add(ndviNorm.multiply(0.35))
  .add(cropNorm.multiply(0.25))
  .rename('Drought_Risk')
  .clip(studyArea);

// ============================================================
// 7️⃣ CLASSIFY into 5 risk levels
// ============================================================
var riskClass = ee.Image(0)
  .where(droughtRisk.gt(0.00).and(droughtRisk.lte(0.20)), 1)  // Very Low
  .where(droughtRisk.gt(0.20).and(droughtRisk.lte(0.40)), 2)  // Low
  .where(droughtRisk.gt(0.40).and(droughtRisk.lte(0.60)), 3)  // Moderate
  .where(droughtRisk.gt(0.60).and(droughtRisk.lte(0.80)), 4)  // High
  .where(droughtRisk.gt(0.80), 5)                              // Very High
  .clip(studyArea)
  .rename('Risk_Class');

// ============================================================
// 8️⃣ VISUALIZE
// ============================================================
var riskVis = {
  min: 0, max: 1,
  palette: ['1a9850','91cf60','ffffbf','fc8d59','d73027']
};

var classVis = {
  min: 1, max: 5,
  palette: ['1a9850','91cf60','ffffbf','fc8d59','d73027']
};

// Individual indicators
Map.addLayer(spiNorm,   {min:0,max:1, palette:['ffffcc','fd8d3c','800026']},
  'SPI Drought Frequency (normalized)', false);
Map.addLayer(ndviNorm,  {min:0,max:1, palette:['006837','ffffcc','d73027']},
  'NDVI Stress (normalized)', false);
Map.addLayer(cropNorm,  {min:0,max:1, palette:['f7fcf5','74c476','00441b']},
  'Crop Exposure (normalized)', false);

// Composite risk
Map.addLayer(droughtRisk, riskVis, 'Drought Risk Index — Philippines 1991–2020');
Map.addLayer(riskClass,   classVis,'Drought Risk Classes (5-level)', false);

// ============================================================
// 9️⃣ LEGEND
// ============================================================
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '10px 14px',
    backgroundColor: 'rgba(255,255,255,0.95)'
  }
});

legend.add(ui.Label({
  value: '🌵 Drought Risk Index | Philippines | 1991–2020',
  style: {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}
}));
legend.add(ui.Label({
  value: 'SPI 40% + NDVI Anomaly 35% + Crop Exposure 25%',
  style: {fontSize: '10px', color: '#666', margin: '0 0 6px 0'}
}));

// Continuous color bar
var colorBar = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select(0)
    .divide(100)
    .visualize({min: 0, max: 1, palette: riskVis.palette}),
  params: {bbox: '0,0,1,0.1', dimensions: '200x20'},
  style: {stretch: 'horizontal', margin: '0 0 4px 0', maxHeight: '20px'}
});
legend.add(colorBar);

var barLabels = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style:{margin:'0 0 8px 0'}});
barLabels.add(ui.Label('Low risk',  {fontSize:'10px', color:'#555'}));
barLabels.add(ui.Label('',          {stretch:'horizontal'}));
barLabels.add(ui.Label('High risk', {fontSize:'10px', color:'#555'}));
legend.add(barLabels);

// Class labels
var riskClasses = [
  {color: '1a9850', label: 'Very Low   (0.00–0.20)'},
  {color: '91cf60', label: 'Low        (0.20–0.40)'},
  {color: 'ffffbf', label: 'Moderate   (0.40–0.60)'},
  {color: 'fc8d59', label: 'High       (0.60–0.80)'},
  {color: 'd73027', label: 'Very High  (0.80–1.00)'}
];

legend.add(ui.Label('Risk classes:', {fontSize:'10px', color:'#666', margin:'2px 0 4px 0'}));

riskClasses.forEach(function(c) {
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
  value: 'Toggle individual indicators in Layers panel',
  style: {fontSize:'9px', color:'#aaa', margin:'6px 0 0 0'}
}));

Map.add(legend);

// ============================================================
// 🔟 ZONAL STATISTICS
// ============================================================

// Area per risk class (km²)
var areaImage = ee.Image.pixelArea().divide(1e6).addBands(riskClass);

var areaStats = areaImage.reduceRegion({
  reducer: ee.Reducer.sum().group({groupField:1, groupName:'class'}),
  geometry: studyArea,
  scale: 5556,
  maxPixels: 1e13,
  bestEffort: true
});
print('📊 Area per drought risk class (km²):', areaStats);

// National mean risk score
print('📈 National mean drought risk score (0–1):', droughtRisk.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: studyArea,
  scale: 5556,
  bestEffort: true
}));

// Mean risk per indicator
print('🌧️ Mean SPI drought frequency (0–1):', spiNorm.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: studyArea,
  scale: 5566,
  bestEffort: true
}));

print('🌿 Mean NDVI stress (0–1):', ndviNorm.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: studyArea,
  scale: 5556,
  bestEffort: true
}));

print('🌾 Mean crop exposure (0–1):', cropNorm.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: studyArea,
  scale: 10,
  bestEffort: true
}));

// ============================================================
// 1️⃣1️⃣ EXPORT TO DRIVE
// ============================================================

// Composite risk index (continuous 0–1)
Export.image.toDrive({
  image: droughtRisk,
  description: 'DroughtRiskIndex_Philippines_1991_2020',
  folder: 'GEE_Exports',
  fileNamePrefix: 'DroughtRiskIndex_Philippines_1991_2020',
  region: studyArea,
  scale: 5556,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Classified risk (5 classes, integer)
Export.image.toDrive({
  image: riskClass,
  description: 'DroughtRiskClass_Philippines_1991_2020',
  folder: 'GEE_Exports',
  fileNamePrefix: 'DroughtRiskClass_Philippines_1991_2020',
  region: studyArea,
  scale: 5556,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// Individual indicators (for validation/reporting)
Export.image.toDrive({
  image: spiNorm.addBands(ndviNorm).addBands(cropNorm),
  description: 'DroughtRisk_Indicators_Philippines_1991_2020',
  folder: 'GEE_Exports',
  fileNamePrefix: 'DroughtRisk_Indicators_Philippines_1991_2020',
  region: studyArea,
  scale: 5556,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
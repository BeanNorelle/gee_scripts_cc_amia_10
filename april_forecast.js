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
// 2️⃣ CFSv2 — TEST with April 2024
// ==============================
var cfsv2 = ee.ImageCollection('NOAA/CFSV2/FOR6H')
  .select('Precipitation_rate_surface_6_Hour_Average')
  .filterDate('2024-04-01', '2024-05-01');

// Check if collection has data using evaluate
cfsv2.size().evaluate(function(size) {
  var hasData = size > 0;
  print('📡 CFSv2 collection has data:', hasData);
  print('📡 CFSv2 collection size:', size);
  
  // Create forecast image with fallback
  var aprilForecastImg;
  if (hasData) {
    aprilForecastImg = cfsv2
      .map(function(img) {
        return img.multiply(21600) // convert to mm per 6-hour step
                  .rename('precip_mm')
                  .copyProperties(img, ['system:time_start']);
      })
      .sum()
      .clip(studyArea)
      .rename('Forecast_mm');
  } else {
    aprilForecastImg = ee.Image.constant(0)
      .clip(studyArea)
      .rename('Forecast_mm');
  }
  
  // ==============================
  // 3️⃣ CHIRPS baseline (1991–2020 April)
  // ==============================
  var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
    .select('precipitation');
  
  function monthlySum(year, month) {
    var start = ee.Date.fromYMD(year, month, 1);
    var end = start.advance(1, 'month');
    return chirps
      .filterDate(start, end)
      .sum()
      .set('year', year);
  }
  
  var years = ee.List.sequence(1991, 2020);
  
  var aprilBaseline = ee.ImageCollection(
    years.map(function(y) {
      return monthlySum(ee.Number(y).toInt(), 4);
    })
  );
  
  var baselineMean = aprilBaseline.mean().clip(studyArea);
  
  // ==============================
  // 4️⃣ Match CHIRPS to CFSv2 grid
  // ==============================
  
  // Get projection - use CFSv2 if available, otherwise CHIRPS
  var targetProj;
  if (hasData) {
    targetProj = ee.Image(cfsv2.first()).select(0).projection();
  } else {
    targetProj = baselineMean.projection();
  }
  
  // Resample CHIRPS to target resolution
  var baselineCoarse = baselineMean
    .resample('bilinear')
    .reproject({
      crs: targetProj,
      scale: 5000
    })
    .rename('Baseline_mm');
  
  // ==============================
  // 5️⃣ Anomaly calculations
  // ==============================
  
  var anomalyImg;
  var anomalyPctImg;
  var droughtClass;
  
  if (hasData) {
    // Calculate anomaly only if we have forecast data
    anomalyImg = aprilForecastImg.subtract(baselineCoarse).rename('Anomaly_mm');
    
    anomalyPctImg = anomalyImg
      .divide(baselineCoarse.where(baselineCoarse.eq(0), 1))
      .multiply(100)
      .rename('Anomaly_pct');
    
    // Drought classification
    droughtClass = anomalyPctImg.expression(
      "(p <= -50) ? 1" +   // severe drought
      ": (p <= -25) ? 2" + // moderate drought
      ": (p < 25) ? 3" +   // near normal
      ": (p < 50) ? 4" +   // moderately wet
      ": 5",               // very wet
      {p: anomalyPctImg}
    ).rename('Drought_Class');
  } else {
    // Create placeholder images when no forecast data
    anomalyImg = ee.Image.constant(0).clip(studyArea).rename('Anomaly_mm');
    anomalyPctImg = ee.Image.constant(0).clip(studyArea).rename('Anomaly_pct');
    droughtClass = ee.Image.constant(3).clip(studyArea).rename('Drought_Class');
  }
  
  // ==============================
  // 6️⃣ Visualization
  // ==============================
  
  var forecastVis = {
    min: 0,
    max: 400,
    palette: [
      'ffffff','c6e9f7','74c8e8',
      '2196f3','1565c0','0d47a1',
      '1b5e20','388e3c','f9a825',
      'e65100','b71c1c'
    ]
  };
  
  var anomalyPctVis = {
    min: -60,
    max: 60,
    palette: [
      'a50026','f46d43','fdae61',
      'ffffff',
      'a6d96a','1a9850','006837'
    ]
  };
  
  var droughtVis = {
    min: 1,
    max: 5,
    palette: [
      '8b0000', // severe drought
      'e31a1c', // moderate drought
      'ffffbf', // normal
      '41b6c4', // wet
      '081d58'  // very wet
    ]
  };
  
  Map.addLayer(aprilForecastImg, forecastVis, 'Forecast Rainfall (mm)', hasData);
  Map.addLayer(anomalyPctImg, anomalyPctVis, 'Anomaly (%)', false);
  Map.addLayer(droughtClass, droughtVis, 'Drought Classification', false);
  
  // ==============================
  // 7️⃣ Regional statistics
  // ==============================
  var regions = ee.FeatureCollection("FAO/GAUL/2015/level1")
    .filter(ee.Filter.eq('ADM0_NAME', 'Philippines'));
  
  if (hasData) {
    var regionalStats = anomalyPctImg.reduceRegions({
      collection: regions,
      reducer: ee.Reducer.mean(),
      scale: 5000
    });
    print('📊 Regional anomaly %:', regionalStats);
    
    // ==============================
    // 8️⃣ Area per drought class
    // ==============================
    var classArea = ee.Image.pixelArea().divide(1e6) // km²
      .addBands(droughtClass)
      .reduceRegion({
        reducer: ee.Reducer.sum().group({
          groupField: 1,
          groupName: 'class'
        }),
        geometry: studyArea,
        scale: 5000,
        maxPixels: 1e13
      });
    print('📊 Area per drought class (km²):', classArea);
    
    // ==============================
    // 9️⃣ Exports
    // ==============================
    Export.image.toDrive({
      image: aprilForecastImg,
      description: 'Forecast_Rainfall_April2024_PH',
      region: studyArea,
      scale: 5000,
      maxPixels: 1e13,
      folder: 'EE_exports'
    });
    
    Export.image.toDrive({
      image: anomalyPctImg,
      description: 'Rainfall_Anomaly_Percent_April2024_PH',
      region: studyArea,
      scale: 5000,
      maxPixels: 1e13,
      folder: 'EE_exports'
    });
    
    Export.image.toDrive({
      image: droughtClass,
      description: 'Drought_Classification_April2024_PH',
      region: studyArea,
      scale: 5000,
      maxPixels: 1e13,
      folder: 'EE_exports'
    });
    
    Export.table.toDrive({
      collection: regionalStats,
      description: 'Regional_Drought_Stats_April2024_PH',
      fileFormat: 'CSV',
      folder: 'EE_exports'
    });
  } else {
    print('⚠️ No forecast data available for April 2024');
    print('💡 Try changing the date to a more recent month');
    print('💡 CFSv2 forecasts are typically available for the current year only');
  }
});
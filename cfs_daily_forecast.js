var cfs = ee.ImageCollection("NOAA/CFSV2/FOR6H")
  .select('Precipitation_rate_surface_6_Hour_Average');

  // rainfall (mm) = rate × seconds
  // 6 hrs = 21600 seconds

  var cfsRain = cfs.map(function(img) {
  return img.multiply(21600) // convert rate → mm per 6 hrs
    .copyProperties(img, ['system:time_start']);
});


function dailySum(date) {
  var start = ee.Date(date);
  var end = start.advance(1, 'day');

  return cfsRain
    .filterDate(start, end)
    .sum()
    .set('system:time_start', start.millis());
}

// Generate daily collection
var startDate = ee.Date('2011-01-01');
var endDate   = ee.Date('2020-12-31');

var nDays = endDate.difference(startDate, 'day');

var dailyCollection = ee.ImageCollection(
  ee.List.sequence(0, nDays.subtract(1)).map(function(d) {
    return dailySum(startDate.advance(d, 'day'));
  })
);


//April rainfall per year

function aprilCFS(year) {
  var start = ee.Date.fromYMD(year, 4, 1);
  var end   = start.advance(1, 'month');

  return dailyCollection
    .filterDate(start, end)
    .sum()
    .set('year', year)
    .set('system:time_start', start.millis());
}


//build hindcast collection

var cfsAprilHindcast = ee.ImageCollection(
  ee.List.sequence(2011, 2020).map(function(y) {
    return aprilCFS(y);
  })
);


print('CFS April Hindcasts:', cfsAprilHindcast);

Map.addLayer(
  cfsAprilHindcast.mean(),
  {min: 0, max: 400, palette: ['white','blue','green','yellow','red']},
  'CFS April Mean'
);
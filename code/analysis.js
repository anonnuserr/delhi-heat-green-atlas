// =====================================================================
// DELHI HEAT–GREEN ATLAS  ·  Analysis script (version 4)
// Green cover, surface heat and greening priority in Delhi
// Author: Ekorshi Paul Choudhury
//
// HOW TO RUN: paste into https://code.earthengine.google.com, click Run.
// The full run takes 1–3 minutes. Results appear in the Console tab,
// maps appear below, exports are in the Tasks tab.
//
// What this script does
//   1. Builds a pre-monsoon (15 Apr – 15 Jun) cloud-free Landsat 8/9
//      composite for 2023, 2024 and 2025.
//   2. Computes vegetation (NDVI) and land surface temperature (LST, °C).
//   3. Keeps only tree cover, grassland and built-up land (ESA WorldCover),
//      so farmland, bare soil and water do not distort the comparison.
//   4. Measures, for each year: the NDVI–LST correlation, °C change per
//      +0.1 NDVI, and mean LST of trees vs grass vs built-up.
//   5. Builds a greening-priority index for built-up land
//      (heat + lack of vegetation + population), maps the top 10 %
//      and lists the 10 largest priority clusters.
//   NEW IN VERSION 4
//   6. Same-neighbourhood test: splits Delhi into 1 km squares and compares
//      trees, grass and built-up land INSIDE each square, so north–south
//      climate differences cannot explain the result.
//   7. Planting opportunity: finds open grassland that is hotter than the
//      average built-up surface and lies within 500 m of a priority
//      hotspot (existing open land where canopy could go, next to the
//      people who need it), with hectares per priority cluster.
// =====================================================================


// ---------- 0. SETTINGS -------------------------------------------------
var YEARS      = [2023, 2024, 2025];
var MAIN_YEAR  = 2025;
var SEASON     = ['-04-15', '-06-15'];   // pre-monsoon window
var MAX_CLOUD  = 20;                     // % scene cloud cover
var N_SAMPLES  = 4000;                   // random sample points per year
var SEED       = 42;


// ---------- 1. STUDY AREA: NCT of Delhi --------------------------------
// FAO GAUL boundary; falls back to a rectangle if the name isn't found.
var rect = ee.Geometry.Rectangle([76.84, 28.40, 77.35, 28.89]);
var gaul = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM0_NAME', 'India'))
  .filter(ee.Filter.inList('ADM1_NAME',
    ['Delhi', 'NCT of Delhi', 'National Capital Territory of Delhi']));
var delhi = ee.Geometry(ee.Algorithms.If(gaul.size().gt(0), gaul.geometry(), rect));
print('1 · Official Delhi boundary used (should be true):', gaul.size().gt(0));


// ---------- 2. LANDSAT 8/9: NDVI + LST ---------------------------------
function prep(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 0).eq(0)      // fill
    .and(qa.bitwiseAnd(1 << 1).eq(0))          // dilated cloud
    .and(qa.bitwiseAnd(1 << 2).eq(0))          // cirrus
    .and(qa.bitwiseAnd(1 << 3).eq(0))          // cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0));         // cloud shadow

  // Surface reflectance (Collection 2 scale factors)
  var sr = img.select(['SR_B4', 'SR_B5']).multiply(0.0000275).add(-0.2);
  var ndvi = sr.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

  // Surface temperature, Kelvin -> °C. ST_B10 = 0 means no data.
  var st = img.select('ST_B10');
  var lst = st.multiply(0.00341802).add(149.0).subtract(273.15).rename('LST');
  var valid = st.gt(0).and(lst.gt(10)).and(lst.lt(75));   // physical sanity check

  return ee.Image(ndvi.addBands(lst)
    .updateMask(clear).updateMask(valid)
    .copyProperties(img, ['system:time_start']));
}

function landsat(year) {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
    .filterBounds(delhi)
    .filterDate(year + SEASON[0], year + SEASON[1])
    .filter(ee.Filter.lt('CLOUD_COVER', MAX_CLOUD))
    .map(prep);
}


// ---------- 3. LAND COVER: ESA WorldCover 2021 (10 m) -----------------
// 10 tree cover · 30 grassland · 50 built-up
// Excluded: 40 cropland, 60 bare/sparse, 80 water (and other classes).
var wc    = ee.ImageCollection('ESA/WorldCover/v200').first();
var trees = wc.eq(10);
var grass = wc.eq(30);
var built = wc.eq(50);
var keep  = trees.or(grass).or(built);

function composite(year) {
  return landsat(year).median().clip(delhi)
    .updateMask(keep)
    .select(['NDVI', 'LST']);
}


// ---------- 4. PER-YEAR STATISTICS -------------------------------------
function meanLST(img, maskImg) {
  return ee.Number(img.select('LST').updateMask(maskImg).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: delhi, scale: 30,
    maxPixels: 1e10, tileScale: 4
  }).get('LST'));
}

function yearStats(year) {
  var col = landsat(year);
  var img = composite(year);
  var s = img.updateMask(img.select('NDVI').gt(0))   // drop any residual water
    .sample({region: delhi, scale: 30, numPixels: N_SAMPLES,
             seed: SEED, tileScale: 4});
  var corr = s.reduceColumns(ee.Reducer.pearsonsCorrelation(), ['NDVI', 'LST']);
  var fit  = s.reduceColumns(ee.Reducer.linearFit(), ['NDVI', 'LST']);
  var t = meanLST(img, trees);
  var g = meanLST(img, grass);
  var b = meanLST(img, built);
  return ee.Feature(null, {
    year: year,
    scenes: col.size(),
    sample_points: s.size(),
    r_NDVI_LST: ee.Number(corr.get('correlation')),
    degC_per_0_1_NDVI: ee.Number(fit.get('scale')).divide(10),
    LST_trees: t,
    LST_grass: g,
    LST_built: b,
    built_minus_trees: b.subtract(t),
    built_minus_grass: b.subtract(g)
  });
}

var summary = ee.FeatureCollection(YEARS.map(yearStats));

print('2 · RESULTS TABLE (one row per year):',
  ui.Chart.feature.byFeature(summary, 'year', [
    'scenes', 'sample_points', 'r_NDVI_LST', 'degC_per_0_1_NDVI',
    'LST_trees', 'LST_grass', 'LST_built',
    'built_minus_trees', 'built_minus_grass'
  ]).setChartType('Table'));

// Dot comparison: trees vs grass vs built-up, per year
print('3 · Mean surface temperature by land cover:',
  ui.Chart.feature.byFeature(summary, 'year', ['LST_trees', 'LST_grass', 'LST_built'])
    .setChartType('ColumnChart')
    .setOptions({
      title: 'Mean pre-monsoon land surface temperature, Delhi',
      vAxis: {title: '°C'}, hAxis: {title: 'Year', format: '####'},
      colors: ['#1b7837', '#b8a14a', '#c2410c']
    }));


// ---------- 5. MAIN-YEAR SCATTER PLOT -----------------------------------
var comp = composite(MAIN_YEAR);
var samples = comp.updateMask(comp.select('NDVI').gt(0)).sample({
  region: delhi, scale: 30, numPixels: N_SAMPLES, seed: SEED,
  geometries: true, tileScale: 4
});

print('4 · Scatter, ' + MAIN_YEAR + ' (each dot = one 30 m pixel):',
  ui.Chart.feature.byFeature(samples, 'NDVI', 'LST')
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Delhi, pre-monsoon ' + MAIN_YEAR + ': vegetation vs surface temperature',
      hAxis: {title: 'NDVI (greenness)'},
      vAxis: {title: 'Land surface temperature (°C)'},
      pointSize: 2, dataOpacity: 0.4, legend: 'none',
      trendlines: {0: {color: '#111111', lineWidth: 2}}
    }));


// ---------- 6. GREENING PRIORITY INDEX (built-up land only) -----------
// Population per 100 m cell, GHSL P2023A, 2025 epoch
var pop = ee.Image('JRC/GHSL/P2023A/GHS_POP/2025').select('population_count');

// Rescale to 0–1 between the 5th and 95th percentile inside Delhi
function stretch(img) {
  var p = img.reduceRegion({
    reducer: ee.Reducer.percentile([5, 95]), geometry: delhi,
    scale: 100, maxPixels: 1e10, bestEffort: true, tileScale: 4
  });
  var v  = ee.List(p.values());          // keys sort as *_p5, *_p95
  var lo = ee.Number(v.get(0));
  var hi = ee.Number(v.get(1));
  return img.subtract(lo).divide(hi.subtract(lo)).clamp(0, 1);
}

var heat     = stretch(comp.select('LST'));                         // hotter     -> 1
var greenGap = ee.Image(1).subtract(stretch(comp.select('NDVI')));  // less green -> 1
var people   = stretch(pop);                                        // more people-> 1

var priority = ee.Image.cat([heat, greenGap, people])
  .reduce(ee.Reducer.mean()).rename('priority')
  .updateMask(built).clip(delhi);

var p90 = ee.Number(ee.List(priority.reduceRegion({
  reducer: ee.Reducer.percentile([90]), geometry: delhi,
  scale: 30, maxPixels: 1e10, bestEffort: true, tileScale: 4
}).values()).get(0));
print('5 · Priority score that marks the top 10 % (0–1):', p90);

var hotspots = priority.gte(p90).selfMask().toInt().rename('hotspot');

// The 10 largest connected hotspot clusters (look these up on Google Maps)
var clusters = hotspots.reduceToVectors({
  geometry: delhi, scale: 100, geometryType: 'polygon',
  eightConnected: true, maxPixels: 1e10, tileScale: 4
}).map(function (f) {
  var c = f.geometry().centroid(10).coordinates();
  return f.set({
    area_ha: f.geometry().area(10).divide(1e4),
    lat: c.get(1),
    lon: c.get(0)
  });
}).sort('area_ha', false).limit(10);

print('6 · Ten largest priority clusters (paste lat, lon into Google Maps):',
  ui.Chart.feature.byFeature(clusters, 'area_ha', ['lat', 'lon'])
    .setChartType('Table'));


// ---------- 6b. SAME-NEIGHBOURHOOD (PAIRED) COMPARISON ------------------
// Delhi is cut into 1 km squares. Only squares containing at least 20
// pixels (≈1.8 ha) each of tree cover, grassland AND built-up land are
// kept. Within each square we compare the three, then summarise across
// squares. This removes city-scale gradients (e.g. the cooler south).
var GRID = delhi.coveringGrid(ee.Projection('EPSG:32643').atScale(1000));
var MIN_PIX = 20;

function pairedStats(year) {
  var lst = composite(year).select('LST');
  var stack = ee.Image.cat([
    lst.updateMask(trees).rename('T'),
    lst.updateMask(grass).rename('G'),
    lst.updateMask(built).rename('B')
  ]);
  var cells = stack.reduceRegions({
    collection: GRID,
    reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.count(), sharedInputs: true}),
    scale: 30, tileScale: 4
  })
  .filter(ee.Filter.gte('T_count', MIN_PIX))
  .filter(ee.Filter.gte('G_count', MIN_PIX))
  .filter(ee.Filter.gte('B_count', MIN_PIX))
  .map(function (f) {
    var t = ee.Number(f.get('T_mean'));
    var g = ee.Number(f.get('G_mean'));
    var b = ee.Number(f.get('B_mean'));
    return f.set({
      g_minus_t: g.subtract(t),
      b_minus_t: b.subtract(t),
      g_minus_b: g.subtract(b),
      trees_coolest: t.lt(g).and(t.lt(b)),
      grass_hotter_than_built: g.gt(b)
    });
  });
  var med = ee.List(cells.reduceColumns(
    ee.Reducer.median().repeat(3), ['g_minus_t', 'b_minus_t', 'g_minus_b']).get('median'));
  return ee.Feature(null, {
    year: year,
    squares_compared: cells.size(),
    median_grass_minus_trees: med.get(0),
    median_built_minus_trees: med.get(1),
    median_grass_minus_built: med.get(2),
    pct_squares_trees_coolest: ee.Number(cells.aggregate_mean('trees_coolest')).multiply(100),
    pct_squares_grass_hotter_than_built: ee.Number(cells.aggregate_mean('grass_hotter_than_built')).multiply(100)
  });
}

var paired = ee.FeatureCollection(YEARS.map(pairedStats));
print('7 · SAME-NEIGHBOURHOOD TEST (1 km squares, °C differences):',
  ui.Chart.feature.byFeature(paired, 'year', [
    'squares_compared', 'median_grass_minus_trees', 'median_built_minus_trees',
    'median_grass_minus_built', 'pct_squares_trees_coolest',
    'pct_squares_grass_hotter_than_built'
  ]).setChartType('Table'));


// ---------- 6c. PLANTING OPPORTUNITY: hot open grassland near hotspots ---
var builtMean2025 = meanLST(comp, built);
var hotGrass = comp.select('LST').gt(builtMean2025).and(grass).selfMask().rename('hot_grass');
var nearHotspot = hotspots.unmask(0).clip(delhi)
  .focalMax({radius: 500, units: 'meters', kernelType: 'circle'}).gt(0);
var opportunity = hotGrass.updateMask(nearHotspot).rename('opportunity');

function hectares(maskImg) {
  return ee.Number(ee.Image.pixelArea().updateMask(maskImg).reduceRegion({
    reducer: ee.Reducer.sum(), geometry: delhi, scale: 30,
    maxPixels: 1e10, tileScale: 4
  }).get('area')).divide(1e4);
}

print('8 · PLANTING OPPORTUNITY, ' + MAIN_YEAR + ':', ee.Dictionary({
  built_up_mean_LST_C: builtMean2025,
  hot_open_grassland_ha_citywide: hectares(hotGrass),
  hot_open_grassland_within_500m_of_hotspots_ha: hectares(opportunity)
}));

var clusterOpp = ee.Image.pixelArea().divide(1e4).updateMask(opportunity).rename('opp_ha')
  .reduceRegions({
    collection: clusters.map(function (f) { return f.buffer(500, 10); }),
    reducer: ee.Reducer.sum(), scale: 30, tileScale: 4
  });
print('9 · Hot open grassland within 500 m of each of the 10 largest clusters (ha):',
  ui.Chart.feature.byFeature(clusterOpp, 'area_ha', ['lat', 'lon', 'sum'])
    .setChartType('Table'));


// ---------- 7. MAPS ------------------------------------------------------
Map.setOptions('HYBRID');
Map.centerObject(delhi, 10);
Map.addLayer(comp.select('NDVI'),
  {min: 0, max: 0.6, palette: ['#e8dcc0', '#a6d96a', '#1a9641', '#00441b']},
  'A · Vegetation (NDVI)', false);
Map.addLayer(comp.select('LST'),
  {min: 35, max: 50, palette: ['#2c7bb6', '#abd9e9', '#ffffbf', '#fdae61', '#d7191c']},
  'B · Surface temperature (°C)', true, 0.85);
Map.addLayer(priority,
  {min: 0.3, max: 0.9, palette: ['#fff7bc', '#fec44f', '#d95f0e', '#7f2704']},
  'C · Greening priority (built-up land)', false);
Map.addLayer(hotspots, {palette: ['#e7298a']}, 'D · Top 10 % hotspots', false);
Map.addLayer(clusters.style({color: '#e7298a', fillColor: '00000000', width: 2}),
  {}, 'E · 10 largest clusters', false);
Map.addLayer(hotGrass, {palette: ['#fdd835']}, 'F · Hot open grassland (hotter than avg built-up)', false);
Map.addLayer(opportunity, {palette: ['#00e5ff']}, 'G · Planting opportunity (hot grass ≤500 m from hotspots)', false);
Map.addLayer(ee.Image().paint(ee.FeatureCollection([ee.Feature(delhi)]), 0, 2),
  {palette: ['#ffffff']}, 'Delhi boundary');


// ---------- 8. EXPORTS (Tasks tab -> Run) -------------------------------
Export.table.toDrive({
  collection: summary, description: 'delhi_results_by_year', fileFormat: 'CSV'
});
Export.table.toDrive({
  collection: clusters, description: 'delhi_priority_clusters', fileFormat: 'CSV'
});
Export.table.toDrive({
  collection: samples, description: 'delhi_samples_' + MAIN_YEAR, fileFormat: 'CSV'
});
Export.table.toDrive({
  collection: paired, description: 'delhi_same_neighbourhood_test', fileFormat: 'CSV'
});
Export.table.toDrive({
  collection: clusterOpp, description: 'delhi_planting_opportunity_by_cluster', fileFormat: 'CSV'
});
Export.image.toDrive({
  image: comp.addBands(priority).toFloat(),
  description: 'delhi_ndvi_lst_priority_' + MAIN_YEAR,
  region: delhi, scale: 30, maxPixels: 1e10
});

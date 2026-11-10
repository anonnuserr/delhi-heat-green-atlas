// =====================================================================
// DELHI HEAT–GREEN ATLAS  ·  Interactive map app
// Author: Ekorshi Paul Choudhury
//
// HOW TO PUBLISH AS A PUBLIC WEBSITE
//   1. Paste this into a NEW script in https://code.earthengine.google.com
//      and click Run to check it works.
//   2. Click "Apps" (top of the editor) -> "New App".
//   3. Name it e.g. delhi-heat-green-atlas, keep "Current contents of
//      editor" as the source, click Publish.
//   4. You get a public link like
//      https://YOUR-PROJECT.projects.earthengine.app/view/delhi-heat-green-atlas
//      (it can take a few minutes to go live).
// =====================================================================

var YEAR = 2025;

// ---------- Study area ---------------------------------------------------
var rect = ee.Geometry.Rectangle([76.84, 28.40, 77.35, 28.89]);
var gaul = ee.FeatureCollection('FAO/GAUL/2015/level1')
  .filter(ee.Filter.eq('ADM0_NAME', 'India'))
  .filter(ee.Filter.inList('ADM1_NAME',
    ['Delhi', 'NCT of Delhi', 'National Capital Territory of Delhi']));
var delhi = ee.Geometry(ee.Algorithms.If(gaul.size().gt(0), gaul.geometry(), rect));

// ---------- Landsat 8/9 NDVI + LST (same method as the analysis) -------
function prep(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 0).eq(0)
    .and(qa.bitwiseAnd(1 << 1).eq(0))
    .and(qa.bitwiseAnd(1 << 2).eq(0))
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0));
  var sr = img.select(['SR_B4', 'SR_B5']).multiply(0.0000275).add(-0.2);
  var ndvi = sr.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
  var st = img.select('ST_B10');
  var lst = st.multiply(0.00341802).add(149.0).subtract(273.15).rename('LST');
  var valid = st.gt(0).and(lst.gt(10)).and(lst.lt(75));
  return ee.Image(ndvi.addBands(lst).updateMask(clear).updateMask(valid)
    .copyProperties(img, ['system:time_start']));
}

var col = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
  .filterBounds(delhi)
  .filterDate(YEAR + '-04-15', YEAR + '-06-15')
  .filter(ee.Filter.lt('CLOUD_COVER', 20))
  .map(prep);

var wc    = ee.ImageCollection('ESA/WorldCover/v200').first();
var trees = wc.eq(10);
var grass = wc.eq(30);
var built = wc.eq(50);
var keep  = trees.or(grass).or(built);

var compAll = col.median().clip(delhi).select(['NDVI', 'LST']);   // every land type
var comp = compAll.updateMask(keep);                               // analysis land types only
var showAll = false;                                               // checkbox state

// ---------- Greening priority (built-up land) ---------------------------
var pop = ee.Image('JRC/GHSL/P2023A/GHS_POP/2025').select('population_count');

function stretch(img) {
  var p = img.reduceRegion({
    reducer: ee.Reducer.percentile([5, 95]), geometry: delhi,
    scale: 100, maxPixels: 1e10, bestEffort: true, tileScale: 4
  });
  var v = ee.List(p.values());
  var lo = ee.Number(v.get(0));
  var hi = ee.Number(v.get(1));
  return img.subtract(lo).divide(hi.subtract(lo)).clamp(0, 1);
}

var priority = ee.Image.cat([
    stretch(comp.select('LST')),
    ee.Image(1).subtract(stretch(comp.select('NDVI'))),
    stretch(pop)
  ]).reduce(ee.Reducer.mean()).rename('priority')
  .updateMask(built).clip(delhi);

var p90 = ee.Number(ee.List(priority.reduceRegion({
  reducer: ee.Reducer.percentile([90]), geometry: delhi,
  scale: 30, maxPixels: 1e10, bestEffort: true, tileScale: 4
}).values()).get(0));
var hotspots = priority.gte(p90).selfMask();

// ---------- Layer definitions -------------------------------------------
var LAYERS = {
  'Surface temperature (°C)': {
    image: comp.select('LST'), imageAll: compAll.select('LST'),
    vis: {min: 35, max: 50, palette: ['#2c7bb6', '#abd9e9', '#ffffbf', '#fdae61', '#d7191c']},
    low: '35 °C', high: '50 °C',
    note: 'Landsat land surface temperature, pre-monsoon ' + YEAR +
          ' (median of clear scenes, 15 Apr – 15 Jun). Surface, not air, temperature.'
  },
  'Vegetation (NDVI)': {
    image: comp.select('NDVI'), imageAll: compAll.select('NDVI'),
    vis: {min: 0, max: 0.6, palette: ['#e8dcc0', '#a6d96a', '#1a9641', '#00441b']},
    low: '0 (bare)', high: '0.6 (dense)',
    note: 'Normalised Difference Vegetation Index from the same Landsat scenes. Higher = greener.'
  },
  'Greening priority': {
    image: priority,
    vis: {min: 0.3, max: 0.9, palette: ['#fff7bc', '#fec44f', '#d95f0e', '#7f2704']},
    low: 'lower', high: 'higher',
    note: 'Built-up land only. Equal-weight average of heat, lack of vegetation and ' +
          'population (each rescaled 0–1 within Delhi).'
  },
  'Top 10 % priority hotspots': {
    image: hotspots,
    vis: {palette: ['#e7298a']},
    low: null, high: null,
    note: 'Built-up pixels in the top 10 % of the greening-priority score.'
  }
};

// ---------- Map ----------------------------------------------------------
Map.setOptions('HYBRID');
Map.centerObject(delhi, 11);
Map.style().set('cursor', 'crosshair');
Map.setControlVisibility({layerList: false, drawingToolsControl: false});

var boundary = ui.Map.Layer(
  ee.Image().paint(ee.FeatureCollection([ee.Feature(delhi)]), 0, 2),
  {palette: ['#ffffff']}, 'Delhi boundary');

function showLayer(name) {
  var L = LAYERS[name];
  var img = (showAll && L.imageAll) ? L.imageAll : L.image;
  allBox.style().set('shown', !!L.imageAll);
  Map.layers().reset([ui.Map.Layer(img, L.vis, name, true, 0.85), boundary]);
  legendPanel.clear();
  legendPanel.add(ui.Label(name, {fontWeight: 'bold', margin: '4px 0 6px 0'}));
  if (L.low !== null) {
    legendPanel.add(ui.Thumbnail({
      image: ee.Image.pixelLonLat().select(0),
      params: {bbox: [0, 0, 1, 0.1], dimensions: '220x12', format: 'png',
               min: 0, max: 1, palette: L.vis.palette},
      style: {stretch: 'horizontal', margin: '0', maxHeight: '18px'}
    }));
    legendPanel.add(ui.Panel([
      ui.Label(L.low,  {margin: '2px 0', fontSize: '11px'}),
      ui.Label('',     {stretch: 'horizontal'}),
      ui.Label(L.high, {margin: '2px 0', fontSize: '11px'})
    ], ui.Panel.Layout.flow('horizontal')));
  } else {
    legendPanel.add(ui.Panel([
      ui.Label('', {backgroundColor: '#e7298a', padding: '7px', margin: '2px 6px 0 0'}),
      ui.Label('Hotspot', {fontSize: '11px', margin: '2px 0'})
    ], ui.Panel.Layout.flow('horizontal')));
  }
  legendPanel.add(ui.Label(L.note, {fontSize: '11px', color: '#555555', margin: '6px 0 0 0'}));
  if (L.imageAll && !showAll) {
    legendPanel.add(ui.Label(
      'Uncoloured areas are cropland, bare land or water. They are left out of the ' +
      'analysis because bare and ploughed ground heats up without being urban. ' +
      'Tick the box above to colour them too.',
      {fontSize: '11px', color: '#555555', margin: '6px 0 0 0'}));
  }
}

// ---------- Side panel ---------------------------------------------------
var panel = ui.Panel({style: {width: '340px', padding: '14px'}});

panel.add(ui.Label('Delhi Heat–Green Atlas',
  {fontSize: '22px', fontWeight: 'bold', margin: '0 0 4px 0'}));
panel.add(ui.Label('Where would trees cool Delhi the most?',
  {fontSize: '14px', margin: '0 0 10px 0'}));
panel.add(ui.Label(
  'Satellite evidence on how green cover relates to surface heat across Delhi, ' +
  'and which built-up neighbourhoods are hot, lack vegetation and are densely ' +
  'populated: candidates for nature-based cooling.',
  {fontSize: '12px', color: '#444444'}));

panel.add(ui.Label('Map layer', {fontWeight: 'bold', margin: '12px 0 2px 0'}));
var selector = ui.Select({
  items: Object.keys(LAYERS),
  value: 'Surface temperature (°C)',
  onChange: showLayer,
  style: {stretch: 'horizontal'}
});
panel.add(selector);

var allBox = ui.Checkbox({
  label: 'Also colour cropland, bare land and water',
  value: false,
  onChange: function (v) { showAll = v; showLayer(selector.getValue()); },
  style: {fontSize: '12px', margin: '6px 0 0 8px'}
});
panel.add(allBox);

var legendPanel = ui.Panel({style: {margin: '8px 0 0 0'}});
panel.add(legendPanel);

panel.add(ui.Label('Click anywhere on the map',
  {fontWeight: 'bold', margin: '14px 0 2px 0'}));
var inspector = ui.Panel([ui.Label('Values for the clicked 30 m pixel appear here.',
  {fontSize: '12px', color: '#555555'})]);
panel.add(inspector);

panel.add(ui.Label(
  'Data: Landsat 8/9 Collection 2 (USGS), ESA WorldCover 2021, ' +
  'GHSL population 2025 (JRC), FAO GAUL boundary. ' +
  'Preliminary exploratory analysis by Ekorshi Paul Choudhury, ' +
  'B.Tech Mathematics & Computing, SMVDU.',
  {fontSize: '10px', color: '#777777', margin: '16px 0 0 0'}));

ui.root.insert(0, panel);

// ---------- Click inspector ----------------------------------------------
var LC_NAMES = {10: 'Tree cover', 20: 'Shrubland', 30: 'Grassland', 40: 'Cropland',
                50: 'Built-up', 60: 'Bare / sparse', 70: 'Snow / ice', 80: 'Water',
                90: 'Wetland', 95: 'Mangroves', 100: 'Moss / lichen'};

function fmt(x, digits, unit) {
  return (x === null || x === undefined) ? 'n/a' : x.toFixed(digits) + (unit || '');
}

Map.onClick(function (coords) {
  inspector.clear();
  inspector.add(ui.Label('Loading…', {fontSize: '12px', color: '#555555'}));
  var pt = ee.Geometry.Point(coords.lon, coords.lat);
  var stack = col.median().select(['NDVI', 'LST'])
    .addBands(priority)
    .addBands(wc.rename('landcover'));
  stack.reduceRegion({reducer: ee.Reducer.first(), geometry: pt, scale: 30})
    .set('inDelhi', delhi.contains(pt, 1))
    .evaluate(function (v, err) {
      inspector.clear();
      if (err || !v) {
        inspector.add(ui.Label('Could not read this point. Try another spot.',
          {fontSize: '12px', color: '#b00020'}));
        return;
      }
      if (!v.inDelhi) {
        inspector.add(ui.Label('This point is outside the Delhi study area (white outline). ' +
          'Click inside the outline to read values.', {fontSize: '12px', color: '#555555'}));
        return;
      }
      var lc = LC_NAMES[v.landcover] || 'n/a';
      var inAnalysis = v.landcover === 10 || v.landcover === 30 || v.landcover === 50;
      var rows = [
        ['Location', coords.lat.toFixed(4) + '° N, ' + coords.lon.toFixed(4) + '° E'],
        ['Land cover', lc],
        ['In the analysis', inAnalysis ? 'Yes' : 'No (only trees, grass and built-up are compared)'],
        ['Surface temperature', fmt(v.LST, 1, ' °C')],
        ['Vegetation (NDVI)', fmt(v.NDVI, 2)],
        ['Greening priority', v.priority === null || v.priority === undefined
            ? 'only scored on built-up land' : fmt(v.priority, 2) + ' (0–1)']
      ];
      rows.forEach(function (r) {
        inspector.add(ui.Panel([
          ui.Label(r[0], {fontSize: '12px', color: '#555555', width: '130px', margin: '2px 0'}),
          ui.Label(r[1], {fontSize: '12px', fontWeight: 'bold', margin: '2px 0'})
        ], ui.Panel.Layout.flow('horizontal')));
      });
    });
});

showLayer('Surface temperature (°C)');

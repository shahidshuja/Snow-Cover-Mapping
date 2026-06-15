// ============================================================
//  Snow Cover Mapping — Multi-Method NDSI Comparison
//  Sensor   : Sentinel-2 SR (10–20 m)
//  Methods  : Standard NDSI | Log-transformed NDSI |
//             SWIR-corrected NDSI
//  Baseline : Three-year temporal mean (2020–2022)
//  Author   : Shahid Shuja Shafai  <shahidshafai@gmail.com>
//  Lab      : Himalayan Cryospheric Research Lab,
//             University of Kashmir
// ============================================================

// ── DESCRIPTION ─────────────────────────────────────────────
//
// Accurate snow cover delineation in high-mountain environments
// is complicated by two persistent challenges: (1) spectral
// confusion between fresh snow and glacial lakes, which share
// similar high reflectance in the visible and (2) SWIR band
// saturation over optically deep snow packs.
//
// This script evaluates three NDSI formulations against a
// three-year (2020–2022) multi-temporal mean Sentinel-2 composite
// to establish a temporally stable snow cover baseline:
//
//  Method 1 — Standard NDSI:
//    (Green − SWIR1) / (Green + SWIR1)
//    Threshold > 0.2 (MODIS/Landsat convention).
//    Effective for broad snow mapping but cannot reliably
//    separate fresh snow from glacial lakes in the Himalayas.
//
//  Method 2 — Log-transformed NDSI (sign-inverted):
//    −1 × (log(Green) − log(SWIR1)) / (log(Green) + log(SWIR1))
//    The −1 sign inversion is intentional: by operating in
//    log-reflectance space, the relative contrast between fresh
//    snow (very high SWIR absorption) and glacial lake surfaces
//    (moderate SWIR reflectance) is amplified. The inversion
//    preserves the convention that high positive values = snow.
//    Scaled ×1000; empirical threshold > 400.
//
//  Method 3 — SWIR-corrected NDSI:
//    Standard NDSI ÷ √(1 / SWIR1)
//    Dividing by √(1/B11) penalises pixels with very low SWIR
//    reflectance (i.e., SWIR-saturated deep snow), compressing
//    the dynamic range in a physically motivated way and reducing
//    false positives over water bodies.
//    Scaled ×1000; empirical threshold > 230.
//
// The three-year mean composite suppresses transient cloud
// artefacts, seasonal noise, and ephemeral snow events, giving
// a robust representation of persistent snow cover across the
// study area. Thresholds (0.2 / 400 / 230) were empirically
// determined from the dataset distribution.
//
// ── HOW TO USE ───────────────────────────────────────────────
//  1. Replace AOI below with your study area geometry or asset.
//  2. Adjust the date window and cloud threshold if needed.
//  3. Run the script — all three methods display simultaneously.
//  4. Compare the binary snow masks to assess which method best
//     separates snow from glacial lakes in your region.
//  5. Use the Tasks tab to export any of the three products.
//
// ─────────────────────────────────────────────────────────────


// ============================================================
//  ★  USER-CONFIGURABLE PARAMETERS — EDIT HERE  ★
// ============================================================

// ── AOI ──────────────────────────────────────────────────────
// Replace with your own study area. Options:
//   a) Drawn geometry:   var AOI = geometry;
//   b) Uploaded asset:   var AOI = ee.FeatureCollection('users/you/boundary');
//   c) Admin boundary:   var AOI = ee.FeatureCollection('FAO/GAUL/2015/level2')
//                                    .filter(ee.Filter.eq('ADM2_NAME','YourDistrict'));
var AOI = geometry;   // <-- REPLACE with your study area

// ── TEMPORAL BASELINE ────────────────────────────────────────
// Three-year window builds a stable multi-annual mean composite.
// Extend or narrow as needed; ensure winter months are included.
var BASELINE_START = '2020-01-01';
var BASELINE_END   = '2022-12-30';

// ── CLOUD COVER THRESHOLD ────────────────────────────────────
// Maximum scene-level cloud cover (%). Lower = fewer but cleaner images.
var MAX_CLOUD_PCT = 20;

// ── NDSI THRESHOLDS ──────────────────────────────────────────
// Empirically determined from the 2020–2022 dataset distribution.
// Adjust if applying to a different region or time period.
var THRESHOLD_STANDARD  = 0.2;    // Standard NDSI (unitless, −1 to 1)
var THRESHOLD_LOG       = 400;    // Log-NDSI       (×1000 scaled)
var THRESHOLD_SWIR      = 230;    // SWIR-corrected (×1000 scaled)

// ── EXPORT SETTINGS ──────────────────────────────────────────
var EXPORT_FOLDER = 'GEE_Outputs';
var EXPORT_SCALE  = 20;   // metres — matches S2 SWIR band resolution

// ============================================================
//  END OF USER PARAMETERS
// ============================================================


// ── Map Initialisation ───────────────────────────────────────
Map.centerObject(AOI, 10);
Map.setOptions('SATELLITE');


// ── Cloud Masking — Sentinel-2 SR ────────────────────────────
// Masks clouds (QA60 bit 10) and cirrus (bit 11), scales DN → [0,1]
function maskS2clouds(image) {
  var qa            = image.select('QA60');
  var cloudBitMask  = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var clearMask = qa.bitwiseAnd(cloudBitMask).eq(0)
                    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(clearMask).divide(10000);
}


// ── Build Multi-Year Mean Composite ─────────────────────────
print('── Building ' + BASELINE_START.slice(0,4) +
      '–' + BASELINE_END.slice(0,4) + ' Sentinel-2 composite ──');

var s2Collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate(BASELINE_START, BASELINE_END)
  .filterBounds(AOI)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_CLOUD_PCT))
  .map(maskS2clouds)
  .select(['B2', 'B3', 'B4', 'B8', 'B9', 'B11', 'B12']);

print('Collection size (cloud-filtered):', s2Collection.size());

var img = s2Collection.mean().clip(AOI);

// True colour reference layer
Map.addLayer(
  img,
  { min: 0.05, max: 0.35, bands: ['B4', 'B3', 'B2'] },
  'True Colour Composite'
);
Map.addLayer(
  img,
  { min: 0.05, max: 0.76, bands: ['B8', 'B4', 'B3'] },
  'NIR False Colour Composite'
);


// ── Method 1: Standard NDSI ──────────────────────────────────
// NDSI = (Green − SWIR1) / (Green + SWIR1)
// Widely used for snow detection (Hall et al. 1995).
// Threshold > 0.2 follows MODIS/Landsat convention.
// Limitation: cannot reliably separate fresh snow from glacial
// lake surfaces in Himalayan terrain.

var ndsiStandard = img.normalizedDifference(['B3', 'B11']).rename('NDSI_Standard');

var snowMaskStandard = ndsiStandard.updateMask(ndsiStandard.gt(THRESHOLD_STANDARD)).toInt();

Map.addLayer(ndsiStandard,    { min: -0.3, max: 0.8 }, 'Method 1 — Standard NDSI (continuous)');
Map.addLayer(snowMaskStandard, {},                      'Method 1 — Standard NDSI Snow Mask (threshold > ' + THRESHOLD_STANDARD + ')');


// ── Method 2: Log-transformed NDSI (sign-inverted) ───────────
// Operates in log-reflectance space to amplify contrast between
// fresh snow (near-zero SWIR reflectance → very negative log) and
// glacial lakes (moderate SWIR → less negative log).
//
// Formula: −1 × (log(B3) − log(B11)) / (log(B3) + log(B11))
//
// The −1 sign inversion is intentional: in log-reflectance space
// the numerator sign flips relative to standard NDSI, so inverting
// restores the convention that high values = snow, while also
// providing enhanced separation of snow from water bodies.
// Scaled ×1000 to improve readability. Threshold > 400.

var logB3  = img.select('B3').log();
var logB11 = img.select('B11').log();
var logImg = logB3.addBands(logB11);

var ndsiLog = logImg.expression(
  '-1 * ((b3 - b11) / (b3 + b11))',
  { b3: logImg.select('B3'), b11: logImg.select('B11') }
).multiply(1000).clip(AOI).rename('NDSI_Log');

var snowMaskLog = ndsiLog.updateMask(ndsiLog.gt(THRESHOLD_LOG)).toInt();

Map.addLayer(ndsiLog,    { min: -221, max: 961 }, 'Method 2 — Log NDSI (continuous, ×1000)');
Map.addLayer(snowMaskLog, {},                      'Method 2 — Log NDSI Snow Mask (threshold > ' + THRESHOLD_LOG + ')');


// ── Method 3: SWIR-corrected NDSI ───────────────────────────
// Divides standard NDSI by √(1/B11) — equivalently multiplies by
// √(B11). This penalises pixels with very low SWIR reflectance
// (SWIR-saturated deep snow), compressing the index range in a
// physically motivated way and reducing false positives over
// spectrally similar water bodies.
//
// Formula: NDSI / √(1/B11)  =  NDSI × √(B11)
// Scaled ×1000. Threshold > 230.

var swirCorrection = ee.Image.constant(1).divide(img.select('B11')).sqrt();

var ndsiSwir = ndsiStandard.expression(
  'ndsi / corr',
  { ndsi: ndsiStandard.select('NDSI_Standard'), corr: swirCorrection }
).multiply(1000).clip(AOI).rename('NDSI_SWIR_Corrected');

var snowMaskSwir = ndsiSwir.updateMask(ndsiSwir.gt(THRESHOLD_SWIR)).toInt();

Map.addLayer(ndsiSwir,    { min: -172, max: 283 }, 'Method 3 — SWIR-corrected NDSI (continuous, ×1000)');
Map.addLayer(snowMaskSwir, {},                      'Method 3 — SWIR-corrected NDSI Snow Mask (threshold > ' + THRESHOLD_SWIR + ')');


// ── Snow Cover Area Statistics ───────────────────────────────
print('── Snow Cover Area Statistics ──');

function snowAreaKm2(snowMask, label) {
  var area = ee.Image.pixelArea()
    .updateMask(snowMask)
    .reduceRegion({
      reducer : ee.Reducer.sum(),
      geometry: AOI,
      scale   : EXPORT_SCALE,
      maxPixels: 1e13,
      bestEffort: true
    });
  var areaKm2 = ee.Number(area.get('area')).divide(1e6).round();
  print(label + ' — Snow-covered area (km²):', areaKm2);
}

snowAreaKm2(snowMaskStandard, 'Method 1 (Standard NDSI)');
snowAreaKm2(snowMaskLog,      'Method 2 (Log NDSI)');
snowAreaKm2(snowMaskSwir,     'Method 3 (SWIR-corrected NDSI)');


// ── Export ───────────────────────────────────────────────────
// Export all three snow masks as separate GeoTIFF files.
// Go to Tasks tab (top right) and click Run for each.

Export.image.toDrive({
  image         : ndsiStandard,
  description   : 'NDSI_Standard',
  fileNamePrefix: 'NDSI_Standard_' + BASELINE_START.slice(0,4) + '_' + BASELINE_END.slice(0,4),
  folder        : EXPORT_FOLDER,
  region        : AOI,
  scale         : EXPORT_SCALE,
  maxPixels     : 1e13,
  fileFormat    : 'GeoTIFF',
  formatOptions : { cloudOptimized: true }
});

Export.image.toDrive({
  image         : ndsiLog,
  description   : 'NDSI_Log_Transformed',
  fileNamePrefix: 'NDSI_Log_' + BASELINE_START.slice(0,4) + '_' + BASELINE_END.slice(0,4),
  folder        : EXPORT_FOLDER,
  region        : AOI,
  scale         : EXPORT_SCALE,
  maxPixels     : 1e13,
  fileFormat    : 'GeoTIFF',
  formatOptions : { cloudOptimized: true }
});

Export.image.toDrive({
  image         : ndsiSwir,
  description   : 'NDSI_SWIR_Corrected',
  fileNamePrefix: 'NDSI_SWIR_' + BASELINE_START.slice(0,4) + '_' + BASELINE_END.slice(0,4),
  folder        : EXPORT_FOLDER,
  region        : AOI,
  scale         : EXPORT_SCALE,
  maxPixels     : 1e13,
  fileFormat    : 'GeoTIFF',
  formatOptions : { cloudOptimized: true }
});

print('── Analysis complete — check Tasks tab to export ──');

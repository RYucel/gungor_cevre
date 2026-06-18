// ==============================================================================
// GÜNGÖR ÇÖPLÜĞÜ — GELİŞMİŞ ÇEVRESEL ETKİ ANALİZİ (v2.0)
// Çoklu-Spektral + Zamansal + Termal + Hava Kalitesi + Hasar Kantifikasyonu
// ==============================================================================

// 1. HEDEF LOKASYON VE ETKİ HALKALARI
var exactLon = 33.42811172485354;
var exactLat = 35.250573118833614;
var aoiPoint = ee.Geometry.Point([exactLon, exactLat]);

// Mesafe bazlı etki gradyanı için halka bölgeler
var coreZone   = aoiPoint.buffer(300);                                        // 0-300m (çekirdek)
var innerZone  = aoiPoint.buffer(700).difference(aoiPoint.buffer(300));       // 300-700m
var outerZone  = aoiPoint.buffer(1500).difference(aoiPoint.buffer(700));      // 700-1500m
var farZone    = aoiPoint.buffer(3000).difference(aoiPoint.buffer(1500));     // 1.5-3km

// Kontrol (referans) bölgesi — benzer arazi, yanmadan uzak ~6 km KB
var controlPoint = ee.Geometry.Point([exactLon - 0.07, exactLat + 0.05]);
var controlZone  = controlPoint.buffer(300);

var aoi = aoiPoint.buffer(3000);

Map.centerObject(aoiPoint, 15);
Map.setOptions('SATELLITE');

// 2. TARİH ARALIĞI (geçmişe dönük gerçek veri)
var startDate = '2023-06-01';
var endDate   = '2025-06-19';

// 3. SENTINEL-2 + GELİŞMİŞ BULUT MASKELEME (QA60 + SCL)
function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBit  = 1 << 10;
  var cirrusBit = 1 << 11;
  var qaMask = qa.bitwiseAnd(cloudBit).eq(0).and(qa.bitwiseAnd(cirrusBit).eq(0));

  // SCL: 3=bulut gölgesi, 8=bulut orta, 9=bulut yüksek, 10=sirrus, 11=kar/buz
  var scl = image.select('SCL');
  var sclMask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));

  return image.updateMask(qaMask.and(sclMask))
              .divide(10000)
              .copyProperties(image, ['system:time_start']);
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40))
  .map(maskS2clouds);

// 4. SPEKTRAL İNDİS PAKETİ
function addIndices(image) {
  var ndvi = image.normalizedDifference(['B8',  'B4']).rename('NDVI'); // Bitki sağlığı
  var nbr  = image.normalizedDifference(['B8',  'B12']).rename('NBR');  // Yanma şiddeti (altın standart)
  var ndbi = image.normalizedDifference(['B11', 'B8']).rename('NDBI');  // Çıplak/insan yapısı
  var ndmi = image.normalizedDifference(['B8',  'B11']).rename('NDMI'); // Nem
  var ndwi = image.normalizedDifference(['B3',  'B8']).rename('NDWI');  // Su içeriği

  // BAI — Yanma Alanı İndeksi (kömürleşmiş toprağa çok duyarlı)
  var red = image.select('B4');
  var nir = image.select('B8');
  var bai = ee.Image(1).divide(
    red.subtract(0.1).pow(2).add(nir.subtract(0.06).pow(2))
  ).rename('BAI');

  return image.addBands([ndvi, nbr, ndbi, ndmi, ndwi, bai]);
}

var s2i = s2.map(addIndices);

// 5. MEDYAN KOMPOZİT VE KATMANLAR
var composite = s2i.median().clip(aoi);

Map.addLayer(composite, {bands: ['B4','B3','B2'], min: 0, max: 0.3, gamma: 1.2},
             '1. Gerçek Renk Kompozit', true);
Map.addLayer(composite.select('NDVI'),
             {min: -0.2, max: 0.8, palette: ['6b0000','ff0000','ffff00','00ff00','003b00']},
             '2. NDVI — Bitki Sağlığı', false);
Map.addLayer(composite.select('NBR'),
             {min: -0.5, max: 0.5, palette: ['1a0000','cc0000','ffcc00','99cc00']},
             '3. NBR — Yanma İzleri', false);
Map.addLayer(composite.select('BAI'),
             {min: 0, max: 300, palette: ['ffffff','fff200','ff8800','ff0000','4a0000']},
             '4. BAI — Yanma Alanı İndeksi', false);
Map.addLayer(composite.select('NDBI'),
             {min: -0.3, max: 0.3, palette: ['0066ff','ffffff','ff0000']},
             '5. NDBI — Çıplak/İnsan Yapısı', false);
Map.addLayer(composite.select('NDMI'),
             {min: -0.2, max: 0.5, palette: ['8b4513','ffe5a0','0066cc']},
             '6. NDMI — Nem', false);

// Halka bölgeleri ve kontrol noktası
Map.addLayer(coreZone,      {color: 'red'},    'Çekirdek (0-300m)', false);
Map.addLayer(innerZone,     {color: 'orange'}, 'İç (300-700m)',     false);
Map.addLayer(outerZone,     {color: 'yellow'}, 'Dış (700-1500m)',   false);
Map.addLayer(controlPoint,  {color: 'lime'},   'Kontrol Noktası',   true);

// 6. ZAMANSAL GRAFİKLER — ÇOK KANITLI
var zones = ee.FeatureCollection([
  ee.Feature(coreZone,  {zone: 'Cekirdek (0-300m)'}),
  ee.Feature(innerZone, {zone: 'Ic (300-700m)'}),
  ee.Feature(outerZone, {zone: 'Dis (700-1500m)'}),
  ee.Feature(farZone,   {zone: 'Uzak (1.5-3km)'})
]);

// 6.1 NDVI zaman serisi (etki gradyanı)
var chartNDVI = ui.Chart.image.seriesByRegion({
  imageCollection: s2i.select('NDVI'),
  regions: zones,
  reducer: ee.Reducer.mean(),
  scale: 10,
  seriesProperty: 'zone',
  xProperty: 'system:time_start'
}).setOptions({
  title: '📊 NDVI Etki Gradyanı (Mesafe ↓ → Tahribat ↓ beklenir)',
  vAxis: {title: 'Ortalama NDVI'},
  hAxis: {title: 'Tarih'},
  lineWidth: 2, pointSize: 3,
  colors: ['#b71c1c', '#ef6c00', '#fdd835', '#2e7d32']
});
print(chartNDVI);

// 6.2 NBR zaman serisi (yanma kanıtı)
var chartNBR = ui.Chart.image.seriesByRegion({
  imageCollection: s2i.select('NBR'),
  regions: zones,
  reducer: ee.Reducer.mean(),
  scale: 10,
  seriesProperty: 'zone',
  xProperty: 'system:time_start'
}).setOptions({
  title: '🔥 NBR — Yanma Şiddeti (Düştükçe yanma artar)',
  vAxis: {title: 'NBR'},
  hAxis: {title: 'Tarih'},
  lineWidth: 2, pointSize: 3,
  colors: ['#b71c1c', '#ef6c00', '#fdd835', '#2e7d32']
});
print(chartNBR);

// 6.3 BAI zaman serisi (kömürleşmiş toprak)
var chartBAI = ui.Chart.image.series({
  imageCollection: s2i.select('BAI'),
  region: coreZone,
  reducer: ee.Reducer.mean(),
  scale: 10,
  xProperty: 'system:time_start'
}).setOptions({
  title: '🔥 BAI — Çekirdek Bölge Yanma Alanı İndeksi',
  vAxis: {title: 'BAI (Yükseldikçe yanma artar)'},
  hAxis: {title: 'Tarih'},
  lineWidth: 2, pointSize: 4,
  colors: ['#7b1fa2']
});
print(chartBAI);

// 6.4 Hedef vs. Kontrol karşılaştırması
var cmpZones = ee.FeatureCollection([
  ee.Feature(coreZone,    {zone: 'Yanma Alani'}),
  ee.Feature(controlZone, {zone: 'Kontrol Alani'})
]);
var chartCmp = ui.Chart.image.seriesByRegion({
  imageCollection: s2i.select('NDVI'),
  regions: cmpZones,
  reducer: ee.Reducer.mean(),
  scale: 10,
  seriesProperty: 'zone',
  xProperty: 'system:time_start'
}).setOptions({
  title: '⚖️ Yanma vs. Kontrol Alanı NDVI',
  vAxis: {title: 'NDVI'}, hAxis: {title: 'Tarih'},
  lineWidth: 2, pointSize: 4,
  colors: ['#b71c1c', '#2e7d32']
});
print(chartCmp);

// 7. AKTİF YANGIN TESPİTİ — FIRMS (MODIS + VIIRS)
var firms = ee.ImageCollection('FIRMS')
  .filterBounds(aoi)
  .filterDate(startDate, endDate);

print('🔥 Toplam FIRMS aktif yangın gözlemi:', firms.size());

// FIRMS tarihlerini listele
var fireDates = firms.aggregate_array('system:time_start')
  .map(function(t){ return ee.Date(t).format('YYYY-MM-dd'); });
print('🔥 Yangın gözlem tarihleri:', fireDates);

var fireFreq = firms.sum().clip(aoi);
Map.addLayer(fireFreq.select('T21'),
             {min: 1, max: 30, palette: ['fff200','ff8800','ff0000','4a0000']},
             '7. FIRMS Aktif Yangın Sıklığı', false);

// 8. SENTINEL-5P HAVA KALİTESİ (NO₂, CO, SO₂, Aerosol-İndeks)
function meanS5P(col, band) {
  return ee.ImageCollection(col).filterBounds(aoi).filterDate(startDate, endDate)
           .select(band).mean().clip(aoi);
}

var no2   = meanS5P('COPERNICUS/S5P/NRTI/L3_NO2',    'tropospheric_NO2_column_number_density');
var co    = meanS5P('COPERNICUS/S5P/NRTI/L3_CO',     'CO_column_number_density');
var so2   = meanS5P('COPERNICUS/S5P/NRTI/L3_SO2',    'SO2_column_number_density');
var aerAI = meanS5P('COPERNICUS/S5P/NRTI/L3_AER_AI', 'absorbing_aerosol_index');

Map.addLayer(no2,   {min: 0, max: 0.0003, palette: ['0066cc','ffff00','ff0000']}, '8. NO₂ Troposferik', false);
Map.addLayer(co,    {min: 0, max: 0.05,   palette: ['0066cc','ffff00','ff0000']}, '9. CO', false);
Map.addLayer(so2,   {min: 0, max: 0.0005, palette: ['0066cc','ffff00','ff0000']}, '10. SO₂', false);
Map.addLayer(aerAI, {min: 0, max: 3,      palette: ['ffffff','fff200','ff8800','4a0000']},
             '11. Aerosol İndeksi (Duman)', false);

// NO₂ zaman serisi (bölgesel hava kalitesi trendi)
var chartNO2 = ui.Chart.image.series({
  imageCollection: ee.ImageCollection('COPERNICUS/S5P/NRTI/L3_NO2')
                    .filterBounds(aoi).filterDate(startDate, endDate)
                    .select('tropospheric_NO2_column_number_density'),
  region: coreZone,
  reducer: ee.Reducer.mean(),
  scale: 7000,
  xProperty: 'system:time_start'
}).setOptions({
  title: '☁️ Çekirdek Bölge NO₂ Zaman Serisi',
  vAxis: {title: 'NO₂ (mol/m²)'}, hAxis: {title: 'Tarih'},
  lineWidth: 2, pointSize: 3, colors: ['#6a1b9a']
});
print(chartNO2);

// 9. LANDSAT-9 YÜZEY SICAKLIĞI (Termal — sıcak nokta kanıtı)
var landsat = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
  .filterBounds(aoi).filterDate(startDate, endDate)
  .map(function(img){
    var lst = img.select('ST_B10').multiply(0.00341802).add(149.0).subtract(273.15).rename('LST_C');
    return lst.copyProperties(img, ['system:time_start']);
  });

var lstMedian = landsat.median().clip(aoi);
Map.addLayer(lstMedian,
             {min: 20, max: 60, palette: ['0000cc','00cc00','ffff00','ff0000','4a0000']},
             '12. Yüzey Sıcaklığı °C (Landsat-9)', false);

var chartLST = ui.Chart.image.series({
  imageCollection: landsat,
  region: coreZone,
  reducer: ee.Reducer.mean(),
  scale: 100,
  xProperty: 'system:time_start'
}).setOptions({
  title: '🌡️ Çekirdek Bölge Yüzey Sıcaklığı',
  vAxis: {title: '°C'}, hAxis: {title: 'Tarih'},
  lineWidth: 2, pointSize: 4, colors: ['#c62828']
});
print(chartLST);

// 10. YANGIN ÖNCESİ / SONRASI HASAR KANTİFİKASYONU
// Bu tarihleri FIRMS gözlemlerine göre ayarlayabilirsiniz
var preStart  = '2023-06-01'; var preEnd  = '2023-12-31';
var postStart = '2024-06-01'; var postEnd = '2024-12-31';

var preImg  = s2i.filterDate(preStart,  preEnd).median();
var postImg = s2i.filterDate(postStart, postEnd).median();

var dNDVI = postImg.select('NDVI').subtract(preImg.select('NDVI')).clip(aoi);
var dNBR  = postImg.select('NBR').subtract(preImg.select('NBR')).clip(aoi);

Map.addLayer(dNDVI, {min: -0.5, max: 0.5, palette: ['b71c1c','ffffff','2e7d32']},
             '13. ΔNDVI (Sonra-Önce)', false);
Map.addLayer(dNBR,  {min: -0.5, max: 0.5, palette: ['b71c1c','ffffff','1565c0']},
             '14. ΔNBR — Yanma Hasarı', false);

// Hasarlı alan hesabı (ΔNDVI < -0.2 ve ΔNBR < -0.1 → orta-yüksek şiddetli hasar)
var damageMask = dNDVI.lt(-0.2).and(dNBR.lt(-0.1)).rename('damage');
var damageArea_m2 = damageMask.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: aoi, scale: 20, maxPixels: 1e10
});

var damageStats = damageMask.addBands(ee.Image.pixelArea().rename('area')).reduceRegion({
  reducer: ee.Reducer.sum().group({groupField: 0}),
  geometry: aoi, scale: 20, maxPixels: 1e10
});

print('⚠️ Tahmini toplam hasarlı alan (m²):', damageArea_m2.get('damage'));

// 11. HASAR HARITASINI EXPORT (Drive)
Export.image.toDrive({
  image: dNBR,
  description: 'Gungor_DNBR_DamageMap',
  folder: 'Gungor_Analizi',
  region: aoi, scale: 20, maxPixels: 1e10,
  crs: 'EPSG:32636'
});

Export.image.toDrive({
  image: composite.select(['NDVI','NBR','BAI','NDBI']),
  description: 'Gungor_Composite_Indices',
  folder: 'Gungor_Analizi',
  region: aoi, scale: 10, maxPixels: 1e10,
  crs: 'EPSG:32636'
});

// 12. KONSOL ÖZET RAPORU
print('════════════════════════════════════════════');
print('  GÜNGÖR ÇÖPLÜĞÜ — ÇEVRESEL ETKİ ÖZETİ');
print('════════════════════════════════════════════');
print('📅 Analiz Dönemi    :', startDate, '→', endDate);
print('🛰️  Sentinel-2 sahnesi:', s2.size());
print('🛰️  Landsat-9 sahnesi :', landsat.size());
print('🔥 FIRMS yangın gözlemi:', firms.size());
print('📊 İndeksler         : NDVI, NBR, BAI, NDBI, NDMI, NDWI');
print('🌫️  Hava Kalitesi     : NO₂, CO, SO₂, Aerosol-İndeks');
print('🌍 Etki halkaları    : 0-300 / 300-700 / 700-1500 / 1500-3000 m');
print('⚖️  Kontrol alanı     : Aktif (Kuzey-Batı)');
print('════════════════════════════════════════════');

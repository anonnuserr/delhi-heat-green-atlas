# Delhi Heat–Green Atlas

**Which kind of green actually cools Delhi, which dense neighbourhoods need it most, and where is there open land to plant canopy next to them?**

A satellite study of green cover, land surface temperature and greening priority in Delhi, using Landsat 8/9 imagery from three pre-monsoon seasons (15 April – 15 June, 2023–2025) in Google Earth Engine.

- **Website:** (https://anonnuserr.github.io/delhi-heat-green-atlas/)
- **Interactive map:** https://elite-rider-296513.projects.earthengine.app/view/delhi-heat-green-atlas

## Key findings

| | 2023 | 2024 | 2025 |
|---|---|---|---|
| Correlation, NDVI vs surface temperature (r) | −0.610 | −0.520 | −0.555 |
| °C change per +0.1 NDVI | −0.94 | −1.09 | −0.86 |
| 1 km squares where tree cover was coolest (of 764) | 89% | 77% | 86% |
| Open grassland minus tree cover, same square (median °C) | +1.41 | +1.60 | +1.27 |
| Built-up minus tree cover, same square (median °C) | +1.22 | +0.87 | +1.05 |

- **Canopy cools, open grass does not.** In the dry pre-monsoon months, open grassland ran hotter than tree cover and slightly hotter than built-up land, both city-wide and within the same neighbourhoods.
- **Priority** (heat + lack of vegetation + population, built-up land only) concentrates in dense settlements in west, north-west, north-east and south Delhi and the old city.
- **Room to plant is scarce where it is needed.** Of about 3,750 ha of open grassland hotter than the average built-up surface (2025), about 710 ha lies within 500 m of a top-10 % priority hotspot.

## Method

1. Landsat 8/9 Collection 2 Level 2 scenes over the NCT of Delhi (FAO GAUL boundary), < 20 % cloud; QA-band masking of fill, cloud, cirrus and shadow.
2. NDVI = (NIR − Red) / (NIR + Red); LST (°C) = ST_B10 × 0.00341802 + 149.0 − 273.15; seasonal median per pixel.
3. ESA WorldCover 2021 used to keep only tree cover, grassland and built-up land.
4. Pearson correlation and linear fit on ~2,500 random pixels per season; mean LST per land-cover class.
5. Same-neighbourhood test: 1 km squares with ≥ 1.8 ha of each class (764 squares); medians of within-square differences.
6. Greening priority = mean of heat, (1 − NDVI) and GHSL 2025 population, each rescaled 0–1 (5th–95th percentile); top 10 % = hotspots.
7. Planting opportunity = grassland hotter than mean built-up LST within 500 m of a hotspot.

## Reproduce

Paste [`code/analysis.js`](code/analysis.js) into the [Earth Engine Code Editor](https://code.earthengine.google.com) and click **Run** (1–5 minutes). Results print in the Console; maps appear below; CSV exports are in the Tasks tab. [`code/map_app.js`](code/map_app.js) is the source of the interactive map.

## Data

USGS Landsat 8/9 Collection 2 Level 2 · ESA WorldCover 10 m 2021 v200 (Zanaga et al., 2022) · GHS-POP R2023A (Schiavina et al., 2023, European Commission JRC) · FAO GAUL 2015 level 1.

## Limitations

Satellite surface temperature (late morning), not air temperature; association, not proof of cause; spatially autocorrelated pixels; equal weights in the priority score; WorldCover "grassland" includes dry open ground and scrub; not all open land is plantable.

## Author

Ekorshi Paul Choudhury, B.Tech Mathematics & Computing, Shri Mata Vaishno Devi University. Preliminary exploratory analysis, November 2026. 

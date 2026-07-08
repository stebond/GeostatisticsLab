# Tutorial: Walker Lake Dataset

This tutorial walks through a complete geostatistical workflow using the
**Walker Lake dataset** (Isaaks and Srivastava, 1989): a 470-point sample of
simulated values (variable **V**, in ppm) from the Walker Lake region of Nevada.
This is the standard teaching dataset for demonstrating variography, kriging and
cross-validation.

**Estimated time:** 25–35 minutes  
**File:** `examples/walker_lake_samples_.dat`

---

## 1. Load the data in QGIS

1. Go to **Layer → Add Layer → Add Delimited Text Layer…**
2. Browse to `examples/walker_lake_samples_.dat`.
3. Set:
   - **File format:** Custom delimiters → Space
   - **X field:** `x`
   - **Y field:** `y`
   - **Geometry CRS:** leave as default (dataset uses local metric coordinates)
4. Click **Add** then **Close**.

The layer `walker_lake_samples_` appears in the Layers panel with 470 features.

---

## 2. Open GeostatisticsLab

Click the **GeostatisticsLab** icon in the QGIS toolbar.

---

## 3. Statistics tab

### 3.1 Descriptive statistics

1. **Point layer:** select `walker_lake_samples_`
2. **Field:** select `V`
3. Click **Run Analysis**

The histogram shows a moderately right-skewed distribution with a few high-value
outliers. Note the clustering of samples in certain areas of the domain — this makes
declustering important for this dataset.

### 3.2 Cell Declustering

1. Switch to the **Data Declustering** sub-tab.
2. Set:
   - **Auto:** true
   - **cell_size_min:** 5
   - **cell_size_max:** 80
   - **num. steps:** 20
   - **Cell Size:** 15 *(starting value; you will adjust after inspecting the curve)*
3. Click **Decluster**

The diagnostic curve (weighted mean vs cell size) is drawn. Identify the cell size
corresponding to the minimum or plateau of the curve — this is typically around
**23 map units** for this dataset. Optionally set **user defined** and set a conveninet **Cell Size** to that value and click
**Decluster** again if needed.

A new layer `walker_lake_samples__decluster` is added to the project, containing all
original columns plus a `decluster_weight` field. The **declustering info box** below
the plot shows the cell size used, the naive (unweighted) mean, and the declustered
(weighted) mean.

Click **Export declustering plot** to save the diagnostic curve.

### 3.3 Normality check

1. Switch to the **Q-Q plot** sub-tab.
2. **Transformed variable:** leave blank for now (we will check the raw variable first).
3. Click **Calculate**

The left Q-Q plot for V shows a clear deviation from normality (heavy right tail),
confirmed by a Shapiro-Wilk p-value << 0.05. A Normal Score Transform is not
strictly required for Ordinary Kriging, but is recommended if you plan to use
Simple Kriging or to compute probability maps.

---

## 4. Variogram tab

Select the **Variogram** tab.  
Set the active layer to **`walker_lake_samples_decluster`** (the declustered layer)
and the field to **`V`**.

### 4.1 Variogram Map

1. Switch to the **Variogram Map** sub-tab.
2. Keep all defaults (18 directions, **Camp. princ. directions** enabled).
3. Click **Calculate**

The polar map shows a mild anisotropy. The detected principal direction is
approximately **N160° (GSLIB convention)** with an anisotropy ratio of about **2:1**
(major to minor).

Click **Export variogram map** to save the figure.

### 4.2 Experimental variogram

1. Switch back to the **Variogram** sub-tab.
2. Set:
   - **Lag mode:** `fixed_lag`
   - **Number of directions:** 4
   - **Azimuth:** 0
   - **use_decl.weights:** ✓ (enable to weight pairs by the declustering weights)
3. Click **Calculate**

The plot shows four directional variograms.

**Optional:** enable **Show var box plot** and **Show h distribution** to inspect both the
semivariance distribution and the pair-distance distribution per lag bin — particularly
useful here given the clustered sampling geometry.

Click **Export variogram plot** to save the plot.

### 4.3 Automatic model fitting

1. Set **Num of autofit struct:** 1
2. Click **Auto Fit**
3. In the autofit settings, just enable the spherical model.

The optimizer fits a **Spherical** model. Check:
- Fit Error (SSE) is reasonably low.
- The model reaches a sill close to the sample variance of V.
- The major range is approximately 1.5–2× the minor range.

Enable **show theoretical model** to overlay the fitted curve on the experimental points.

---

## 5. Kriging tab

### 5.1 Kriging setup

1. **Kriging type:** Ordinary (default)
2. **Neighbourhood:**
   - Num. Points Min: 4
   - Num. Points Max: 16
   - Range max: set to the major range from the variogram model
   - Range min: set to the minor range
   - Azimuth: 160
3. **Grid:**
   - Default extent from the data bounding box
   - Nx: 50, Ny: 50

### 5.2 Cross-validation

1. Click **Calculate** (Cross kriging section)

Inspect the results:
- **Mean error** close to 0 → no systematic bias.
- **Cross Krigin Plot:** estimated vs observed scatter around the 1:1 line.
- **Cross kriging map:** check for any spatial pattern in errors that might indicate
  an inappropriate neighbourhood or variogram model.
- **Cross-kriging neighbourhood:** use the **Previous/Next** buttons to inspect the
  search ellipse and neighbor weights for individual data points, especially those
  with large estimation errors.

Click **Export cross-kriging plot** and **Export cross-kriging map** to save the figures.

**Optional — neighbourhood optimisation:**  
Click **pushButton_ck_opt** to run the Latin Hypercube Sampling optimisation. With
470 data points and 80 LHS trials, the search takes a few minutes. The optimised
parameters are written back to the neighbourhood widgets automatically.

### 5.3 Comparison: Ordinary vs Universal Kriging

Run the estimation twice to compare:

**Ordinary Kriging (OK):**
- Output name: `V_ok.tif`
- Click **Calculate!**

**Universal Kriging (UK):**
- Select **Universal** in the Kriging type group.
- Output name: `V_uk.tif`
- Click **Calculate!**

Compare the two estimate rasters visually in QGIS. UK accounts for the linear
spatial trend present in V (increasing values toward the NE corner of the domain)
and typically produces smoother estimates in trend-affected areas.

Compare the cross-validation SSE for OK and UK to decide which model is more
appropriate for this dataset.

### 5.4 Grid estimation output

Each estimation produces a two-band GeoTIFF:
- **Band 1:** kriging estimate (ppm)
- **Band 2:** kriging variance (ppm²)

---

## 6. Probability Maps tab

Compute the probability of V exceeding **500 ppm**, a hypothetical high-value
threshold for this demonstration.

1. **Estimated data:** `V_ppm_ok.tif`
2. **Variance data:**  `V_ppm_var.ok.tif`
3. **threshold value:** 500
4. **direction:** higher
5. Click **Run!**

The probability raster highlights areas where extremely high V values are most likely
— broadly coinciding with the high-value cluster visible in the original data location
map.

---

## 7. Expected outputs

| File | Content |
|---|---|
| `walker_lake_samples__decluster` | In-memory layer: original data + declustering weights |
| `V_ok.tif` and V_ok_var.tif | Ordinary Kriging estimate + variance (ppm) |
| Probability raster | P(V > 500 ppm) per pixel |

---

## Key observations from this dataset

- **Declustering matters:** the naive mean of V overestimates the true spatial average
  because of preferential sampling in high-value areas. The declustered mean is lower.
- **Anisotropy is mild (~2:1):** a three-structure model is needed to capture the
  gradual rise and levelling of the variogram.
- **UK vs OK:** Universal Kriging reduces the estimation variance in areas far from
  the data when a linear trend is present, but may overfit if the trend is not
  well-supported by the data density.

---

## Reference

Isaaks, E.H., Srivastava, R.M. (1989). *An Introduction to Applied Geostatistics*.
Oxford University Press, New York, 561 pp.

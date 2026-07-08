# Tutorial: SIC97 Rainfall Dataset

This tutorial walks through a complete geostatistical workflow using the **SIC97 dataset**
(Dubois, 1998): 100 daily rainfall measurements collected in Switzerland on 8 May 1986,
in units of 1/10th of mm.

**Estimated time:** 20–30 minutes  
**File:** `examples/sic97.dat`

---

## 1. Load the data in QGIS

1. In QGIS, go to **Layer → Add Layer → Add Delimited Text Layer…**
2. Browse to `examples/sic97.dat.csv`.
3. Set:
   - **File format:** Custom delimiters → Space
   - **X field:** `x`
   - **Y field:** `y`
   - **Geometry CRS:** leave as default (the dataset uses Swiss metric coordinates)
4. Click **Add** then **Close**.

The layer `sic97` appears in the Layers panel as a point layer with 100 features.

---

## 2. Open GeostatisticsLab

Click the **GeostatisticsLab** icon in the QGIS toolbar.  
The plugin dialog opens on the **Statistics** tab.

---

## 3. Statistics tab

### 3.1 Descriptive statistics

1. **Point layer:** select `sic97`
2. **Field:** select `value`
3. Click **Run Analysis**

Expected results:

| Statistic | Value |
|---|---|
| Count | 467 |
| Mean | ~184 (1/10 mm) |
| Std dev | ~112 (1/10 mm) |
| Min | 0 |
| Max | 585 |

The histogram shows a moderately right-skewed distribution.  
Click **Export histogram** to save the plot if needed.

### 3.2 Normal Score Transform

The right skew makes a Normal Score Transform advisable before kriging.

1. Leave **New vector layer** blank (the name `value_NSC_TRSF` is suggested automatically).
2. Click **Transform**.

A new in-memory layer `sic97_NSC_TRSF` appears in the Layers panel. Its rows are sorted
in ascending order of the original `value` column — this sorted table will be used
automatically by the back-transform step later.

### 3.3 Check normality (Q-Q plot)

1. Switch to the **Q-Q plot** sub-tab.
2. **Transformed variable:** select `value_NSC_TRSF`
3. Click **Calculate**

The right panel should show a Q-Q plot for the transformed variable that is nearly
linear, with a Shapiro-Wilk p-value > 0.05, confirming satisfactory normality.  
The left panel shows the original variable for comparison.

---

## 4. Variogram tab

Select the **Variogram** tab. Make sure the active layer is `sic97_NSC_TRSF` and the
field is `value_NSC_TRSF`.

### 4.1 Variogram Map

Before computing directional variograms, use the variogram map to identify
the principal anisotropy direction.

1. Switch to the **Variogram Map** sub-tab.
2. Keep all defaults (18 directions, **Camp. princ. directions** enabled).
3. Click **Calculate**.

The polar map shows an elongated pattern oriented roughly **NW–SE**, indicating
anisotropy. The arrows mark the principal directions automatically detected via the
moment-of-inertia method — note the major direction, approximately **N40° (GSLIB
convention)**.

Click **Export variogram map** to save the figure.

### 4.2 Experimental variogram
option A
1. Switch back to the **Variogram** sub-tab.
2. Set:
   - **Lag mode:** `adaptive_lag`
   - **Number of directions:** 2
   - **Azimuth:** 40 (to align with the identified anisotropy)
3. Click **Calculate**

The plugin computes lag spacing automatically using the AutoVar `datspac_horz`
algorithm. Expected values: lag ≈ 8560 m.

The plot shows two directional variograms. The NW–SE direction (major) reaches
the sill at a larger range than the perpendicular direction (minor), confirming
geometric anisotropy with a ratio of approximately **5:1**.
Option B
2. Set:
   - **Lag mode:** `fixed_lag`
   - **Number of directions:** 4
   - **Azimuth:** 0
3. Click **Calculate**

**Optional:** enable **Show box plot** and **Show h distribution** to visualise the
dispersion of both semivariance values and pair distances per lag bin.

Click **Export variogram** to save the plot.

### 4.3 Automatic model fitting

1. Set **Num of autofit struct:** 2
2. Click **Auto Fit**

The optimizer (L-BFGS-B with multi-start azimuth) fits a two-structure **NUGGET+CUBIC**
model. The **Fit Error (SSE)** label updates with the weighted residual.

Check that:
- The model curve passes through the experimental points reasonably well.
- The sill is close to 1.0 (expected for Normal Score data).
- The range in the major direction is approximately twice that in the minor direction.

Enable **show theoretical model** if it is not already active.

---

## 5. Kriging tab

### 5.1 Kriging setup

1. **Kriging type:** Ordinary
2. **Neighbourhood:**
   - Num. Points Min: 4
   - Num. Points Max: 16
   - Range max: set to the major range from the variogram model
   - Range min: set to the minor range
   - Azimuth: 40
3. **Grid:**
   - Use the default extent (auto-computed from the data bounding box)
   - Nx: 50, Ny: 32 (or adjust to the desired resolution)

### 5.2 Cross-validation

Before running the full grid estimation, validate the model.

1. Click **Calculate** (Cross kriging section)

The cross-validation table fills with observed vs estimated values for all 467 points.
Check the **Summary** row at the bottom:

- **Mean error** should be close to 0 (no systematic bias).
- **SSE** provides a measure of overall prediction error.

Switch to **Cross Krigin Plot** to inspect the scatter of estimated vs observed values;
points should cluster around the 1:1 line.  
Switch to **Cross kriging map** to check for spatial patterns in the errors.

**Optional — neighbourhood optimisation:**  
Click **pushButton_ck_opt** to run the Latin Hypercube Sampling optimisation (80 trials
+ warm-start). The plugin searches for the neighbourhood configuration that minimises
the cross-validation SSE while ensuring all 467 points are estimated.

### 5.3 Grid estimation

1. The output name is suggested automatically (e.g. `value_NSC_TRSF_ok.tif`).
2. Click **Calculate!**

A two-band GeoTIFF is added to the QGIS project:
- **Band 1:** kriging estimate (in Normal Score space)
- **Band 2:** kriging variance

### 5.4 Back-transform

To obtain the estimate in the original rainfall units:

1. **Raster input:** select `value_NSC_TRSF_ok.tif`
2. **Vector layer:** select `sic97_NSC_TRSF`
3. **variable:** `value_NSC_TRSF`
4. Set an output filename (e.g. `sic97_ok_backtransformed.tif`)
5. Click **Back transform**

The back-transformed raster is added to the project in original units (1/10 mm).

---

## 6. Probability Maps tab

Compute the probability of daily rainfall exceeding **200 (1/10 mm)**, a threshold
of practical interest for the Chernobyl deposition study context.

1. **Estimated data:** select Band 1 of `value_NSC_TRSF_ok.tif`
2. **Variance data:** select Band 2 of `value_NSC_TRSF_ok.tif`
3. **threshold value:** 200 *(note: this applies in Normal Score space; convert
   200 (1/10 mm) to its Normal Score equivalent using the transformation table,
   or work directly with the back-transformed raster and use the OK variance as
   an approximation)*
4. **direction:** higher
5. Click **Run!**

The probability map is added to the project. Areas with high probability correspond
to regions where heavy rainfall — and therefore high radioactive deposition — was
most likely.

---

## 7. Expected outputs

| File | Content |
|---|---|
| `sic97_NSC_TRSF` | In-memory layer: Normal Score transformed values + transformation table |
| `value_NSC_TRSF_ok.tif` | Two-band raster: OK estimate + variance (NS space) |
| `sic97_ok_backtransformed.tif` | Estimate in original units (1/10 mm) |
| Probability raster | P(rainfall > 200) per pixel |

---

## Reference

Dubois, G. (1998). Spatial interpolation comparison 97: foreword and introduction.
*Journal of Geographic Information and Decision Analysis* 2(2), 1–10.

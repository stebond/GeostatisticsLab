# GeostatisticsLab

**A open-source QGIS plugin for complete geostatistical workflows**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![QGIS: 3.x](https://img.shields.io/badge/QGIS-3.x-green.svg)](https://www.qgis.org)
[![Python: 3.x](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org)

GeostatisticsLab integrates a **complete geostatistical workflow** into a single QGIS graphical interface, from exploratory data analysis to probability mapping — with no need to manage intermediate files between separate programs.

Developed at **DICAM – University of Bologna**.

---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Plugin Tabs](#plugin-tabs)
  - [Statistics](#statistics-tab)
  - [Variogram](#variogram-tab)
  - [Kriging](#kriging-tab)
  - [Probability Maps](#probability-maps-tab)
- [Example Datasets](#example-datasets)
- [Reproducing the Paper Results](#reproducing-the-paper-results)
- [Repository Structure](#repository-structure)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Features

| Module | Capabilities |
|---|---|
| **Statistics** | Descriptive statistics, histogram, Q-Q plot, Shapiro-Wilk & Anderson-Darling normality tests |
| **Cell Declustering** | GSLIB-style cell declustering with diagnostic curve; declustering weights exported to a new layer |
| **Normal Score Transform** | GSLIB-compatible NST with epsilon jitter; output layer doubles as back-transform table |
| **Variogram** | Directional experimental variogram with **adaptive lag spacing** (replicating AutoVar `datspac_horz`); weighted estimator; variogram map with MOI-based anisotropy detection; 3D surface |
| **Model Fitting** | Manual nested model editor (Nugget, Spherical, Exponential, Gaussian, Linear, Cubic); **Autofit** via L-BFGS-B with multi-start azimuth; GSLIB export |
| **Kriging** | Simple / Ordinary / Universal Kriging on regular grids or point layers; elliptic octant neighborhood search; two-band GeoTIFF output (estimate + variance) |
| **Cross-Validation** | Leave-one-out cross-validation; scatter plot; spatial error map; **interactive neighborhood visualization** (weights + search ellipse per point) |
| **Back-Transform** | Pixel-wise inverse Normal Score Transform using the session transformation table |
| **Probability Maps** | Gaussian exceedance / non-exceedance probability maps from kriging estimate and variance |

---

## Requirements

### Software
- **QGIS 3.16 or later** (tested up to QGIS 3.36)  
  Download: https://www.qgis.org/download

### Python dependencies
All dependencies ship with the standard QGIS Python environment on Windows, macOS and Linux:

| Package | Version | Purpose |
|---|---|---|
| `numpy` | ≥ 1.20 | Numerical arrays and linear algebra |
| `scipy` | ≥ 1.6 | L-BFGS-B optimizer, normality tests |
| `matplotlib` | ≥ 3.3 | All in-plugin plots |
| `osgeo.gdal` | ≥ 3.0 | Raster I/O and in-memory raster creation |

> **Note:** If `scipy` or `matplotlib` are missing from your QGIS Python environment (rare on standard installations), install them via the **OSGeo4W Shell** (Windows) or your system package manager.

### Hardware
- RAM: 4 GB minimum, 8 GB recommended
- Disk: 50 MB for the plugin; additional space for output rasters
- Display: 1280 × 800 or higher

---

## Installation

### Option A — Install from ZIP (recommended)

1. Download the latest release ZIP from the [Releases](../../releases) page.
2. Open QGIS → **Plugins → Manage and Install Plugins → Install from ZIP**.
3. Select the downloaded ZIP file and click **Install Plugin**.
4. The plugin icon will appear in the QGIS toolbar.

### Option B — Manual installation

1. Clone or download this repository:
   ```bash
   git clone https://github.com/[username]/GeostatisticsLab.git
   ```
2. Copy the `GeostatisticsLab/` folder to your QGIS plugins directory:
   - **Windows:** `%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\`
   - **Linux/macOS:** `~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/`
3. Open QGIS → **Plugins → Manage and Install Plugins → Installed** → enable **GeostatisticsLab**.

### Verify installation

After enabling the plugin, its icon should appear in the QGIS toolbar. Click it to open the main dialog.

---

## Quick Start

1. **Load a point layer**: `Layer → Add Layer → Add Delimited Text Layer…` — select one of the example `.dat` files in `examples/`.
2. **Open GeostatisticsLab** by clicking its toolbar icon.
3. In the **Statistics** tab: select the layer and field, press **Run Analysis**.
4. In the **Variogram** tab: select **adaptive_lag**, press **Calculate**, then press **Auto Fit**.
5. In the **Kriging** tab: choose **Ordinary**, set the grid, press **Calculate!**
6. In the **Probability Maps** tab: select the output rasters, set a threshold, press **Run!**

A complete step-by-step tutorial is provided in [`docs/tutorial_sic97.md`](docs/tutorial_sic97.md).

---

## Plugin Tabs

### Statistics Tab

**Inputs**
- Point layer (any QGIS-supported format)
- Numeric field to analyze

**Controls and outputs**

| Control | Description |
|---|---|
| `Run Analysis` | Computes sample count, mean, median, min, max, std dev; draws histogram (15 bins) |
| `Show cumulative histogram` | Overlays empirical CDF on a secondary Y axis |
| `Transformed variable / Calculate` (Q-Q plot sub-tab) | Side-by-side Q-Q plots; Shapiro-Wilk and Anderson-Darling test results |
| `Cell Size`, `cell_size_min/max/steps` / `Decluster` | Runs GSLIB cell declustering; plots weighted-mean vs cell-size diagnostic curve; creates `<layer>_decluster` in-memory layer with `decluster_weight` field |
| `New vector layer` / `Transform` | Applies Normal Score Transform; creates sorted in-memory layer used as back-transform table |

---

### Variogram Tab

**Inputs**
- Point layer (can be the declustered layer from Statistics)
- Numeric field

**Controls and outputs**

| Control | Description |
|---|---|
| `fixed lag` / `adaptive_lag` | Manual or automatic lag spacing (AutoVar `datspac_horz` algorithm) |
| `MaxDistance`, `num_lags`, `lag` | Distance parameters; `lag_tolerance` is always `lag/2` (GSLIB convention) |
| `azimuth`, `Number of directions` | Reference azimuth and number of equally spaced directions over 180° |
| `angular_tolerance`, `bandwidth` | Sector half-aperture and maximum lateral width |
| `use_decl.weights` | Weighted semivariogram using `decluster_weight` column |
| `Show num pairs`, `Show box plot` | Pair counts and within-bin distributions overlaid on plot |
| `show theoretical model` | Overlays current model curve on experimental plot |
| `Normalize by variance` | Displays relative variogram (sill → 1) |
| `Effective distance` | Positions each bin at its actual mean pair distance (AutoVar style) |
| `Calculate` | Computes multi-directional experimental variogram |
| `add / remove structure` | Edits nested theoretical model table (Nugget, Spherical, Exponential, Gaussian, Linear, Cubic) |
| `Auto Fit` | L-BFGS-B optimization with multi-start azimuth; reports Fit Error (SSE) |
| `Settings` sub-tab | Enables/disables model types for Autofit |
| `Export variogram model` | Saves model to GSLIB-format text file |
| `Calculate` (Variogram Map sub-tab) | 2D polar variogram map with MOI-based principal direction arrows |
| `Varigram Model` (3D sub-tab) | 3D semivariance surface with optional theoretical model overlay |

---

### Kriging Tab

**Inputs**
- Point layer with values
- Theoretical variogram model (from Variogram tab)

**Controls and outputs**

| Control | Description |
|---|---|
| `Simple` / `Ordinary` / `Universal` | Kriging type (SK requires known mean; UK includes linear drift) |
| `Num. Points Min/Max` | Minimum and maximum neighbors in the search ellipse |
| `Range max / Range min` | Major and minor semi-axes of the search ellipse |
| `Azimuth` | Orientation of the search ellipse (geographic degrees) |
| `Octant search` / `Min points per octant` | Divides ellipse into 8 sectors with a per-sector minimum |
| `Grid` / `points` radio button | Estimation on a regular grid or on an existing point layer |
| `X/Y Min/Max`, `Nx/Ny`, `dx/dy` | Grid extent and resolution (parameters auto-linked) |
| `Calculate!` | Runs kriging; produces two-band GeoTIFF (estimate + variance) |
| `Calculate` (Cross kriging) | Leave-one-out cross-validation; populates table, scatter plot, error map |
| `<-Previous` / `Next->` | Navigates cross-validation results; shows weights and search ellipse per point |
| `Back transform` | Applies inverse NST to a kriged NS raster using the session transformation table |

**Output raster**
- Band 1: kriging estimate
- Band 2: kriging variance
- NoData value: −9999

---

### Probability Maps Tab

**Inputs**
- Estimated raster (band 1 of kriging output)
- Variance raster (band 2 of kriging output)

**Controls and outputs**

| Control | Description |
|---|---|
| `threshold value` | Threshold for the exceedance/non-exceedance calculation |
| `higher` / `Lower` | Direction: P(Z > threshold) or P(Z ≤ threshold) |
| `Run!` | Computes pixel-wise Gaussian probability map; adds result to QGIS project |

**Formula**
```
P[Z(u) > t] = 1 − Φ{[t − ŷ(u)] / σ_k(u)}    [higher]
P[Z(u) ≤ t] =     Φ{[t − ŷ(u)] / σ_k(u)}    [lower]
```
where Φ is the standard normal CDF and σ_k is the square root of the kriging variance.

---

## Example Datasets

Two public benchmark datasets are included in `examples/`:

### `sic97.dat`
- **Source:** Dubois, G. (1998). Spatial interpolation comparison 97. *Journal of Geographic Information and Decision Analysis* 2(2), 1–10.
- **Content:** 200 measurements of Caesium-137 deposition (Bq/m²) after the Chernobyl accident, southern Germany.
- **Columns:** `x` (Easting, m), `y` (Northing, m), `value` (Bq/m²)
- **Delimiter:** space
- **Recommended field:** `value`

### `walker_lake_samples_.dat`
- **Source:** Isaaks, E.H., Srivastava, R.M. (1989). *An Introduction to Applied Geostatistics*. Oxford University Press.
- **Content:** 470-point sample of the Walker Lake dataset, Nevada.
- **Columns:** `x`, `y`, `V` (primary variable, ppm), `U` (secondary variable)
- **Delimiter:** space
- **Recommended field:** `V`

---

## Reproducing the Paper Results

The results in [*GeostatisticsLab: An open-source QGIS plugin for complete geostatistical workflows* (Computers & Geosciences, in review)] can be reproduced as follows.

### SIC97 case study

1. Load `examples/sic97.dat` as a Delimited Text Layer (X field: `x`, Y field: `y`).
2. **Statistics tab:** field = `value`; press **Run Analysis**.
3. **Statistics → Normal Score Transform:** press **Transform** (output field: `value_NSC_TRSF`).
4. **Variogram tab:** layer = `sic97.dat_NSC_TRSF`; field = `value_NSC_TRSF`; select **adaptive_lag**; Number of directions = 4; press **Calculate**.
5. Press **Auto Fit** (Num of autofit struct = 2).
6. **Kriging tab:** type = Ordinary; Num. Points Min = 4, Max = 16; Range max/min from variogram ranges; Grid: 200 × 200 cells; press **Calculate!**
7. **Probability Maps tab:** select the estimate and variance rasters; threshold = 200; direction = higher; press **Run!**

Expected results: principal direction ≈ N40° (GSLIB); anisotropy ratio ≈ 5:1; two-structure Spherical model.

### Walker Lake case study

1. Load `examples/walker_lake_samples_.dat` (X field: `x`, Y field: `y`).
2. **Statistics → Data Declustering:** cell_size_min = 5, cell_size_max = 80, steps = 20; Cell Size = 20; press **Decluster**.
3. **Variogram tab:** layer = `walker_lake_samples__decluster`; field = `V`; enable **use_decl.weights**; select **adaptive_lag**; press **Calculate**.
4. Press **Auto Fit** (Num of autofit struct = 3).
5. **Kriging tab:** type = Ordinary; Grid: 50 × 50 cells; press **Calculate!**
6. **Cross-Validation:** press **Calculate**; inspect scatter plot and error map.

Expected results: principal direction ≈ N160° (GSLIB); anisotropy ratio ≈ 2:1; three-structure Spherical model.

---

## Repository Structure

```
GeostatisticsLab/
│
├── GeostatisticsLab/                  # QGIS plugin package
│   ├── __init__.py
│   ├── geostatisticslab.py            # Plugin entry point (QGIS Plugin Builder)
│   ├── geostatisticslab_dialog.py     # Main dialog: all algorithms and GUI logic
│   ├── geostatisticslab_dialog_base.ui  # Qt Designer UI file
│   ├── metadata.txt                   # QGIS plugin metadata
│   ├── resources.py                   # Compiled Qt resources
│   └── icons/                         # Plugin icons
│
├── examples/                          # Example datasets
│   ├── sic97.dat                      # SIC97 Caesium-137 dataset (Dubois, 1998)
│   └── walker_lake_samples_.dat       # Walker Lake dataset (Isaaks & Srivastava, 1989)
│
├── docs/                              # Documentation
│   ├── tutorial_sic97.md              # Step-by-step tutorial: SIC97 dataset
│   ├── tutorial_walker_lake.md        # Step-by-step tutorial: Walker Lake dataset
│   └── GeostatisticsLab_UserGuide.docx  # Complete user guide
│
├── tests/                             # Unit tests
│   └── test_algorithms.py             # Tests for core numerical routines
│
├── README.md                          # This file
└── LICENSE                            # GNU General Public License v3
```

---

## Citation

If you use GeostatisticsLab in your research, please cite:

```
[Surname], S. (in review). GeostatisticsLab: An open-source QGIS plugin for complete
geostatistical workflows. Computers & Geosciences.
```

The adaptive lag spacing algorithm implemented in the Variogram tab is described in:

```
Davila Saavedra, L., Deutsch, C.V., 2025. Automatic variogram calculation and modeling.
Computers & Geosciences 195, 105774. https://doi.org/10.1016/j.cageo.2024.105774
```

---

## License

GeostatisticsLab is released under the **GNU General Public License v3.0**.  
See [LICENSE](LICENSE) for the full license text.

---

## Contact

**Stefano [Surname]**  
DICAM – Department of Civil, Chemical, Environmental and Materials Engineering  
University of Bologna, Bologna, Italy  
[email address]

For bug reports and feature requests, please open an [issue](../../issues) on GitHub.

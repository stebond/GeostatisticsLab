# GeostatisticsLab
A complete geostatistical workflow for QGIS: EDA, declustering, Normal Score Transform, variogram analysis and automatic fitting, Simple/Ordinary/Universal Kriging, cross-validation with neighbourhood optimization, back-transform and probability maps. No additional dependencies beyond a standard QGIS installation.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![QGIS: 3.16+](https://img.shields.io/badge/QGIS-3.16+-green.svg)](https://www.qgis.org)
[![Python: 3.x](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org)

Developed at **DICAM – University of Bologna**.

---

## Features

| Tab | What you can do |
|---|---|
| **Statistics** | Descriptive statistics, histogram with cumulative CDF, Q-Q plot with Shapiro-Wilk and Anderson-Darling normality tests, cell declustering with diagnostic curve, Normal Score Transform |
| **Variogram** | Directional experimental variogram with adaptive lag spacing (AutoVar algorithm), horizontal box plots of pair distances per lag bin, manual and automatic model fitting (L-BFGS-B, multi-start azimuth), 2D polar variogram map with anisotropy detection, 3D variogram surface, GSLIB export |
| **Kriging** | Simple / Ordinary / Universal Kriging on regular grids or point layers, elliptic neighbourhood search with optional octant partitioning, leave-one-out cross-validation with scatter plot, spatial error map and interactive neighbourhood inspection, automated neighbourhood parameter optimization via Latin Hypercube Sampling, back-transform of Normal Score rasters |
| **Probability Maps** | Gaussian exceedance / non-exceedance probability maps from kriging estimate and variance |

---

## Requirements

- **QGIS 3.16 or later**
- Python packages already bundled in any standard QGIS installation: `numpy`, `scipy`, `matplotlib`, `osgeo.gdal`

No additional installation is required.

---

## Installation

### From ZIP (recommended)
1. Download the latest release ZIP from the [Releases](../../releases) page.
2. In QGIS: **Plugins → Manage and Install Plugins → Install from ZIP**.
3. Select the ZIP and click **Install Plugin**.

### Manual
1. Clone this repository:
   ```bash
   git clone https://github.com/[username]/GeostatisticsLab.git
   ```
2. Copy the `GeostatisticsLab/` folder to your QGIS plugins directory:
   - **Windows:** `%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\`
   - **Linux/macOS:** `~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/`
3. In QGIS: **Plugins → Manage and Install Plugins → Installed** → enable **GeostatisticsLab**.

---

## Quick Start

1. Load a point layer: **Layer → Add Layer → Add Delimited Text Layer** (use one of the files in `examples/`).
2. Open the plugin from the toolbar.
3. **Statistics tab** — select the layer and field, press **Run Analysis**.
4. **Variogram tab** — select **adaptive_lag**, press **Calculate**, then **Auto Fit**.
5. **Kriging tab** — choose **Ordinary**, configure the grid, press **Calculate!**
6. **Probability Maps tab** — select the output rasters, set a threshold, press **Run!**

Full step-by-step tutorials are in [`docs/tutorial_sic97.md`](docs/tutorial_sic97.md) and [`docs/tutorial_walker_lake.md`](docs/tutorial_walker_lake.md).

---

## Example Datasets

Two public benchmark datasets are included in `examples/`:

| File | Description | Source |
|---|---|---|
| `sic97.dat` | 467 daily rainfall measurements in Switzerland (8 May 1986), 1/10 mm | Dubois (1998) |
| `walker_lake_samples_.dat` | 470-point sample of the Walker Lake dataset, variables V and U | Isaaks & Srivastava (1989) |

### Reproducing the paper results

**SIC97:**
```
1. Load examples/sic97.dat (X: x, Y: y, field: value)
2. Statistics -> Normal Score Transform
3. Variogram -> adaptive_lag, 4 directions -> Calculate -> Auto Fit (2 structures)
4. Kriging -> Ordinary, 50x32 grid -> Calculate!
5. Probability Maps -> threshold 200, direction: higher -> Run!
```
Expected: principal direction ~N40 deg (GSLIB), anisotropy ratio ~5:1, two-structure Spherical model.

**Walker Lake:**
```
1. Load examples/walker_lake_samples.dat (X: x, Y: y, field: V)
2. Statistics -> Decluster (cell size ~20)
3. Variogram -> adaptive_lag, use declustering weights -> Calculate -> Auto Fit (3 structures)
4. Kriging -> Ordinary, 50x50 grid -> Calculate! -> Cross-Validation
```
Expected: principal direction ~N160 deg (GSLIB), anisotropy ratio ~2:1, three-structure Spherical model.

---

## Repository Structure

```
GeostatisticsLab/
├── GeostatisticsLab/
│   ├── __init__.py
│   ├── geostatisticslab.py            # Plugin entry point
│   ├── geostatisticslab_dialog.py     # Main dialog: all algorithms and GUI logic
│   ├── geostatisticslab_dialog_base.ui
│   ├── metadata.txt
│   └── resources.py
├── examples/
│   ├── sic97.dat
│   └── walker_lake_samples_.dat
├── docs/
│   ├── tutorial_sic97.md
│   ├── tutorial_walker_lake.md
│   └── GeostatisticsLab_UserGuide.docx
├── README.md
└── LICENSE
```

---

## Citation

If you use GeostatisticsLab in your research, please cite:

```
Bonduà, S. (in review). GeostatisticsLab: An open-source QGIS plugin for complete
geostatistical workflows. Computers & Geosciences.
```

The adaptive lag spacing algorithm is described in:

```
Davila Saavedra, L., Deutsch, C.V., 2025. Automatic variogram calculation and modeling.
Computers & Geosciences 195, 105774. https://doi.org/10.1016/j.cageo.2024.105774
```

---

## License

Released under the **GNU General Public License v3.0** — see [LICENSE](LICENSE).

---

## Contact

**Stefano Bonduà**
DICAM – University of Bologna
stefano.bondua@unibo.it

Bug reports and feature requests: please open an [issue](../../issues).

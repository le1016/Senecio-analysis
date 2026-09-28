# Senecio-analysis

Analysis code for self-incompatibility studies of *Senecio scandens*.

This repository contains the Python and R scripts used for figure preparation
and expression-data processing for the accompanying Senecio manuscript.
The original scripts are preserved as supplied. Input datasets, manuscript
documents, and figure artwork are not included in this code upload.

## Scripts

| Script | Purpose | Input | Output |
| --- | --- | --- | --- |
| `betabychr.py` | Plot genome-wide beta values in vertically stacked chromosome panels | `Chrall_bt.txt`: tab-separated, no header, columns `CHR START END BETA` | `beta_by_chr.png`, `beta_by_chr.svg` |
| `plot_chr_scatter.py` | Plot raw BetaScan beta values for one chromosome and highlight values above a threshold | The same four-column beta table; chromosome and threshold can be passed as arguments | `<chrom>_betascan_scatter_thr<threshold>.png` and `.pdf` |
| `plot_smoothed_ne.py` | Plot SMC++ effective population-size histories, with optional display smoothing, split times, and interval bands | SMC++ plot CSV or constant-model `model.final.json` files; optional split-time and bootstrap tables | A plot file selected with `--out` (default: `ne_smoothed.pdf`) |
| `Word01.py` | Plot sampling locations over shaded East Asian terrain, with province highlights | `GRAY_HR_SR_W.tif` and `species_points07.csv` with `lon`, `lat`, and optional `species` columns | `Word04.png`, `Word04.pdf`, and an interactive plot |
| `featurecounts_TPM.R` | Convert featureCounts read counts to TPM and calculate average TPM | A featureCounts table with gene lengths and sample counts | `merged.Counts.txt`, `SV_avg_tpm.xls`, `SV_tpm` (tab-separated text) |
| `heatmap.R` | Plot gene-expression heatmaps | `Slocus_tpm03`: a tab-separated table with `Geneid` and four numeric expression columns | `heatmap01.pdf`, `heatmap01.png`, and an interactive heatmap |

## Dependencies

Install the Python packages imported by the scripts:

```bash
python -m pip install -r requirements.txt
```

`scipy` is optional: `plot_smoothed_ne.py` uses it for spline smoothing and
falls back to a moving median if it is unavailable. The geographic plotting
dependencies (`Cartopy`, `rasterio`, and `mplcursors`) are used by `Word01.py`.
That script may download Natural Earth map data through Cartopy when the
required datasets are not already cached.

For the R heatmap script, install:

```r
install.packages(c("pheatmap", "RColorBrewer", "dplyr"))
```

`featurecounts_TPM.R` uses base R. Package versions are not pinned because
the original analysis environment was not supplied with the scripts.

## Running the scripts

Place the required inputs in the working directory, or adjust their paths
in the scripts before running.

### Beta plots

```bash
python betabychr.py
python plot_chr_scatter.py Chrall_bt.txt Chr07 26.44518
```

The genome-wide script sets its threshold in the `THRESHOLD` variable
(currently `33.9388`). The single-chromosome script accepts its threshold
as the third positional argument. These are the values in the supplied
scripts; choose the appropriate analysis threshold for your data.

### Population-size histories

```bash
python plot_smoothed_ne.py --help
python plot_smoothed_ne.py --csv 3pop.marginal.csv --out ne_smoothed.pdf
```

Use `--no-smooth` to disable display smoothing. Spline-model results should
be supplied using the CSV mode; JSON mode requires a model with `epochs`.
CSV time values are used as supplied, while `--gen-time` rescales time in
JSON mode. Check time units before adding split-time lines or glaciation
shading, which are specified in years.

### Sampling map

```bash
python Word01.py
```

Provide a georeferenced DEM covering the configured extent (80–140 degrees
east, 20–45 degrees north) and a sampling-point CSV. The script calls
`plt.show()` after saving the map, so interactive display requires a
graphical environment.

### Expression processing and heatmaps

Before running `featurecounts_TPM.R`, replace its original absolute input
path with the location of your featureCounts table and check the hard-coded
sample-column names against your table. Before running `heatmap.R`, update
its `setwd()` path and check the input filename and expression columns.

```bash
Rscript featurecounts_TPM.R
Rscript heatmap.R
```

The heatmap script reads a separately prepared `Slocus_tpm03` table; it does
not automatically consume the complete TPM output from the other R script.
Its saved heatmaps use the transposed, column-standardized expression data.

## Reproducibility status

All four Python scripts passed a syntax check during repository preparation.
Full figure reproduction has not been verified because the input datasets
are not included in this folder. R scripts were not executed during upload.

## 中文说明

本仓库保存千里光（*Senecio scandens*）自交不亲和研究手稿相关的绘图和
表达数据处理代码，共 4 个 Python 脚本和 2 个 R 脚本。原始代码保持不变。
运行前请准备各脚本所需数据，并修改 R 脚本中的本地路径及样本列名。
手稿、原始数据和排版图件不包含在此次代码上传中。

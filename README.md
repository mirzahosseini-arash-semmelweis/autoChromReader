# autoChromReader

**Automated extraction, quality control, and database-ready organization of retention data from HPLC chromatograms in R**

`autoChromReader` is an R-based workflow for automated processing of large collections of routine chromatograms. It provides a reproducible pipeline from heterogeneous chromatographic files to standardized retention-time and peak-property tables suitable for statistical analysis, QSRR, and machine-learning applications.

The workflow combines vendor-independent chromatogram import, automated baseline correction, statistically calibrated peak detection, optional peak deconvolution, chromatographic calculations, quality-control plots, cross-wavelength consistency checks, and selective operator-guided correction.

## Features

* **Vendor-independent chromatogram import**

  * Uses [`chromConverter`](https://github.com/ethanbass/chromConverter) for supported vendor formats.
  * Directly reads common text exports including CSV, TSV, TXT, DAT, and ARW files.
  * Normalizes chromatograms to a common `time_min` / `intensity` representation.
  * Metadata parsing is modular and can be adapted to local filename conventions.

* **Automated baseline correction**

  * **FlatFit** is the default baseline estimator.
  * FlatFit uses local slope and curvature to identify flat baseline regions and fits a smooth penalized least-squares baseline.
  * The Savitzky--Golay derivative window can be estimated automatically from the characteristic chromatographic feature scale.
  * Optional asymmetric refinement reduces baseline overestimation in densely populated peak regions.
  * **MILE** (morphological and iterative local extremum) is retained as an alternative baseline method.

* **Noise-aware peak detection**

  * Robust noise estimation from relatively quiet chromatographic regions.
  * Topographic prominence-based peak detection.
  * Block-bootstrap calibration of the prominence threshold.
  * User-defined absolute SNR threshold.
  * Optional relative-height threshold to exclude analytically irrelevant minor components.
  * Optional Ricker continuous-wavelet-transform detection of shoulder-like features.

* **Overlapping peak treatment**

  * Automatically identifies clusters of neighboring peaks.
  * Optional exponential-Gaussian hybrid (EGH) deconvolution.
  * Conservative component-count checks prevent deconvolution from silently redefining peak identity.
  * EGH-derived retention time, FWHM, area, asymmetry parameters, and fit diagnostics are retained.

* **Chromatographic calculations**

  * Retention time, $t_R$
  * Dead time, $t_0$
  * Retention factor, $k$
  * FWHM
  * Peak height
  * Peak area
  * SNR
  * Topographic prominence
  * Resolution, $R_s$
  * Selectivity, $\alpha$

* **Quality control**

  * Optional baseline-only inspection pass before peak extraction.
  * Annotated QC plot for every chromatogram.
  * Explicit QC flags for missing peaks, missing $t_0$, saturation, problematic FWHM, and deconvolution failures.
  * Cross-wavelength comparison when the same experiment is recorded at multiple detector wavelengths.
  * Optional cross-wavelength consensus alignment of peak identities.

* **Manual review without repeating the full analysis**

  * Automatically generates an editable `peakfinder_lookup.csv`.
  * Operator can specify:

    * $t_0$
    * approximate `t1`, `t2`, ... positions
    * expected number of peaks
    * exclusions
    * notes
  * Only chromatograms containing operator guidance are reprocessed.

* **Parallel batch processing**

  * Independent chromatograms can be processed using multiple R workers.
  * Designed for datasets containing hundreds of routine chromatographic runs.

---

## Workflow

```text
Chromatographic files
        │
        ▼
read_chrom.R
Vendor-independent import
+ metadata parsing
        │
        ▼
Baseline correction
FlatFit [default] / MILE
        │
        ▼
pick_peaks.R
Noise estimation
+ prominence detection
+ optional CWT shoulders
        │
        ▼
Cluster detection
+ optional EGH deconvolution
        │
        ▼
peakfinder.R
Peak selection
+ chromatographic calculations
+ QC
        │
        ├────► QC plots
        │
        ├────► manual lookup / selective update
        │
        ▼
Database-ready retention table
```

## Repository structure

```text
autoChromReader/
│
├── read_chrom.R
├── baseline_flatfit.R
├── baseline_mile.R
├── pick_peaks.R
├── peakfinder.R
└── README.md
```

### `read_chrom.R`

Handles chromatogram discovery, import, normalization, and metadata parsing.

The standardized data representation contains variables including:

```text
record_uid
logical_key
date
folder_label
filename
file_hash
wavelength_nm
experimentalID
compoundname
columnID
eluentID
modifierID
temp_C
flow_mL_min
time_min
intensity
```

Metadata parsing is intentionally separated from file decoding. Users with different filename conventions can therefore supply their own parser without modifying the chromatographic analysis.

### `baseline_flatfit.R`

Implements FlatFit baseline estimation based on the method introduced in MOCCA2, together with extensions intended for long and densely sampled chromatograms.

Flatness weights are determined from local slope and curvature.

The default implementation estimates the Savitzky--Golay window from the chromatographic feature scale and can optionally perform a small number of asymmetric refinement steps.

### `baseline_mile.R`

Implements morphological and iterative local extremum (MILE) baseline correction as an alternative to FlatFit.

The procedure combines morphological opening, local-extremum selection, PCHIP interpolation, iterative knot refinement, and optional mollifier smoothing.

### `pick_peaks.R`

Contains the main signal-processing functions for:

* robust noise estimation;
* topographic prominence;
* block-bootstrap prominence calibration;
* optional Ricker-CWT shoulder detection;
* peak integration and FWHM estimation;
* peak clustering;
* EGH peak fitting and component selection.

### `peakfinder.R`

Provides the user-facing workflow, including:

```r
inspect_chrom_folder()
analyze_chrom_folder()
update_with_lookup()
make_retention_table()
append_retention_database()
```

---

## Installation

The workflow currently consists of R source files rather than an installed R package.

Required packages can be installed with:

```r
install.packages(c(
  "data.table",
  "stringr",
  "ggplot2",
  "digest",
  "future",
  "future.apply",
  "progressr"
))
```

For vendor-specific chromatographic formats, install `chromConverter` as appropriate.

Then source the workflow:

```r
source("read_chrom.R")
source("baseline_flatfit.R")
source("baseline_mile.R")
source("pick_peaks.R")
source("peakfinder.R")
```

---

## Quick start

A complete folder can be processed with:

```r
result <- analyze_chrom_folder(
  path = "20260814_example",
  baseline_method = "flatfit",
  look_for_t0 = TRUE,
  t0_upper_limit = 1.0,
  snr_min = 20,
  relative_height_min = 0.25,
  n_workers = 16
)
```

FlatFit is the recommended and default baseline method:

```r
baseline_method = "flatfit"
```

MILE can instead be selected with:

```r
baseline_method = "mile"
```

The principal output tables are then available as:

```r
result$runs
result$peaks
result$retention
```

---

## Inspect baselines before peak extraction

For a large new dataset, it can be useful to inspect the imported chromatograms and fitted baselines before running peak detection:

```r
inspect_chrom_folder(
  path = "20260814_example"
)
```

This creates graphical output that can be reviewed for file-reading or baseline-estimation problems.

---

## Output structure

### `result$runs`

One row per chromatogram.

Contains experimental metadata and run-level diagnostics such as:

```text
t0_value_min
t0_source
n_peaks
noise_sd
snr_min_absolute
relative_height_min
effective_analyte_snr_min
prominence_cutoff
baseline_method
baseline_info
qc_status
qc_note
peakfinder_version
```

### `result$peaks`

Long-format peak table containing one row per detected $t_0$ or analyte peak.

Representative variables include:

```text
peak_type
peak_index
peak_source
tR_min
height
prominence
area
SNR
fwhm_min
fwhm_method
saturated
picker_cluster
picker_resolution
deconvolved
deconv_sigma
deconv_tau
k
```

### `result$retention`

Wide, database-oriented representation containing variables such as:

```text
t0

tR1
tR2
...

fwhm1
fwhm2
...

height1
height2
...

area1
area2
...

k1
k2
...

Rs12
Rs23
...

alpha12
alpha23
...
```

Cross-wavelength QC variables are added when applicable.

---

## Cross-wavelength quality control

When chromatograms from the same experimental run are recorded at several detector wavelengths, retention behavior provides an independent internal consistency check.

The workflow compares:

* detected peak counts;
* $t_0$ values;
* corresponding analyte retention times.

For example:

```r
cross_wavelength_tolerance_min = 0.05
```

defines the maximum permitted retention-time spread before the experiment is flagged for review.

An optional consensus-alignment procedure can additionally correct inconsistent **peak identities** across wavelengths. Importantly, this does not average or overwrite measured retention times. Instead, peaks occurring within the specified retention-time tolerance are matched across wavelengths and their `peak_index` assignments are aligned when supported by a majority of the available wavelengths.

This is particularly useful when a weak early peak disappears at one wavelength and would otherwise shift subsequent assignments from `tR1` to `tR2`, etc.

If no sufficiently supported or unambiguous consensus exists, the original assignment is retained and the case is flagged for manual review.

---

## Operator-guided correction

The first automated analysis can create:

```text
peakfinder_lookup.csv
```

The operator can add approximate retention-time hints:

```text
filename,t0,t1,t2,t3,expected_npeaks,exclude,note
```

For example:

```text
sample01.arw,,4.31,5.02,,2,,check weak first peak
```

The reviewed chromatograms can then be selectively reprocessed:

```r
result <- update_with_lookup(result)
```

Only files containing actionable lookup information are recalculated.

---

## Parallel processing

For larger datasets:

```r
result <- analyze_chrom_folder(
  path = "chromatograms",
  n_workers = 16,
  progress = TRUE
)
```

Each worker uses a single `data.table` thread to reduce CPU oversubscription.

In a representative dataset of approximately 400 chromatograms, complete processing—including file import, baseline correction, bootstrap-calibrated peak detection, chromatographic calculations, quality-control assessment, and generation of individual PNG QC plots—required approximately 1 h using 16 workers on an Intel Core Ultra 9-285K workstation with 64 GB DDR5 memory.

Runtime will depend on signal length, bootstrap settings, frequency of overlapping peaks, requested graphical output, and available hardware.

---

## Intended scope

`autoChromReader` is primarily intended for routine analytical chromatographic separations containing a relatively limited number of expected components.

The workflow is particularly suited to:

* chromatographic screening studies;
* retention-time database construction;
* enantioseparation datasets;
* column/mobile-phase screening;
* QSRR datasets;
* machine-learning dataset generation;
* automated extraction of retention parameters from large experimental series.

It is **not intended as a replacement for specialized metabolomics software** or for applications dominated by very complex chromatograms, severe multicomponent coelution, extensive deconvolution, or quantitatively demanding peak-area integration.

For these applications, dedicated tools such as XCMS, MZmine, MOCCA/MOCCA2, or other specialized analytical software may be more appropriate.

---

## Reproducibility and auditability

The workflow deliberately keeps several levels of information separate:

1. original imported chromatograms;
2. processed chromatograms and fitted baselines;
3. automatically detected candidate peaks;
4. final selected peaks;
5. run-level QC information;
6. database-ready retention variables;
7. operator-supplied corrections.

The intention is that automated processing remains inspectable rather than becoming a black-box conversion from raw chromatograms to numerical retention values.

---

## Citation

If you use `autoChromReader` in published work, please cite the associated publication:

> **[Citation to be added upon publication]**

The workflow additionally builds upon methods and software including FlatFit/MOCCA2, MILE, `chromConverter`, continuous-wavelet peak detection, and exponential-Gaussian hybrid peak modeling. Please refer to the manuscript and repository references for the relevant primary publications.

---

## Status

`autoChromReader` is under active development. The current implementation has been developed primarily for automated extraction of retention parameters from routine HPLC chromatograms. Users applying the workflow to substantially different chromatographic data should validate baseline and peak-detection settings for their own application.

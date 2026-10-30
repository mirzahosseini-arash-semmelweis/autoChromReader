# Set working directory to folder containing this R file

source("read_chrom.R")
source("baseline.R")
source("pick_peaks.R")
source("peakfinder.R")

path <- "20260826_sample_data"

# Create raw plots for inspection
inspect_chrom_folder(path)

# Normal first pass that creates lookup table
result <- analyze_chrom_folder_chunked(
  path        = path,
  look_for_t0 = FALSE,
  # set for reproducibility
  max_peaks           = 2L,
  relative_height_min = 0.25,
  do_shoulders        = FALSE,
  min_peak_resolution = 0.05

  # n_workers is set to 1 as default and safeguard ceiling of (n(cores) - 1) is used
  # set n_workers depending on your CPU
)

# ---------------------------------------------------------------------------
# Options (all optional; every resolved value is stored per run in
# result$runs$hyperparameters and in result$settings)
#
# Baseline
#   baseline_method = "flatfit" (default) | "mile" | "psalsa" | "arpls" | "asls"
#   whittaker_args  = list(lambda = NULL, lambda_scale = 4, k = NULL, p = NULL)
#     psalsa/arpls/asls: lambda NULL = automatic from the widest peak
#     (baseline cut-off period = lambda_scale x widest peak base width).
#
# Denoising (after baseline correction, before peak detection)
#   denoise_method  = "none" (default) | "sg" | "gwma" | "ma" | "loess" |
#                     "whittaker" | "spline" | "fft"
#   denoise_args    = list(tune = "acf")   # or "rv", or list(param = 7)
#     The parameter is tuned per chromatogram (Bosten et al. 2023) unless a
#     fixed `param` is given; the noise threshold is re-derived through the
#     same filter, so false-positive control is unchanged.
#
# Shoulders and decomposition
#   do_shoulders              = TRUE   # CWT shoulder search (FALSE = faster)
#   min_peak_resolution       = 0.25   # merge candidates closer than this Rs
#   decomp_require_same_count = TRUE   # FALSE: add BIC-proposed hidden peaks
#   decomp_k_extra            = 2      # hidden components tried per cluster
#   decomp_model_tol          = 0.01   # BEMG model-error floor (x height)
#   decomp_min_fit_rs         = 0.35   # reject fits with components closer
#
# Example for noisy data:
#   result <- analyze_chrom_folder_chunked(path, look_for_t0 = FALSE,
#     max_peaks = 17, baseline_method = "psalsa", denoise_method = "sg")
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Handling t0 marker
# With visible t0 marker use
#    look_for_t0    = TRUE`
#    `t0_upper_limit = 1.5`
# In this case t0 values in lookup act as hints
# For experiments with no t0 marker use
#    `look_for_t0    = FALSE`
# Manually give either one `t0_manual` or per run t0 in lookup.
# Optionally `analyte_lower_limit = 0.8` informs peak picking.
# ---------------------------------------------------------------------------

# After QC and editing lookup only problematic exports are updated
# Run with:
# result <- update_with_lookup(
#   result
# )

# Save results
saveRDS(
  result,
  file.path(path, "peakfinder_result.rds")
)

# Append to database
database <- append_retention_database(
  result$retention,
  database_file = "retention_database.rds"
)

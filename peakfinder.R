## ---------------------------------------------------------------------------
##  peakfinder.R -- Folder-level HPLC retention/FWHM extraction pipeline.
## ---------------------------------------------------------------------------

suppressWarnings(
  suppressPackageStartupMessages({
    library(data.table)
    library(stringr)
    library(ggplot2)
    library(digest)
    library(parallel)
    library(future)
    library(future.apply)
    library(progressr)
  })
)

PEAKFINDER_VERSION <- "0.4.0"

# ---------- Lookup ------------------------------------------------------------

write_peak_lookup_template <- function(
    data,
    path,
    filename  = "peakfinder_lookup.csv",
    max_peaks = 6L,
    overwrite = FALSE
  ) {

  path     <- normalizePath(path, winslash = "/", mustWork = TRUE)
  out_file <- file.path(path, filename)

  if (file.exists(out_file) && !overwrite) {
    message("Lookup file already exists; leaving it unchanged: ", out_file)
    return(invisible(out_file))
  }

  runs <- unique(data[, .(
    record_uid,
    filename,
    date,
    wavelength_nm,
    experimentalID,
    compoundname,
    columnID,
    eluentID,
    modifierID,
    temp_C,
    flow_mL_min
  )])
  
  setorder(runs, filename, wavelength_nm)

  template <- runs[, .(
    record_uid,
    filename,
    wavelength_nm,
    experimentalID,
    compoundname
  )]
  
  setorder(
    template,
    filename,
    wavelength_nm
  )
  
  template[, t0 := NA_real_]
  
  for (j in seq_len(max_peaks)) {
    template[, (paste0("t", j)) := NA_real_]
  }
  
  template[, `:=`(
    expected_npeaks = NA_integer_,
    exclude         = NA,
    note            = NA_character_
  )]

  fwrite(template, out_file, na = "")
  invisible(out_file)
}

read_peak_lookup <- function(path, filename = "peakfinder_lookup.csv") {
  f <- file.path(path, filename)
  if (!file.exists(f)) {
    stop("Lookup file does not exist: ", f)
  }

  x <- fread(f, na.strings = c("", "NA"))
  if (!"record_uid" %in% names(x)) {
    stop("Lookup table must contain a 'record_uid' column.")
  }
  if (anyDuplicated(x$record_uid)) {
    stop("Lookup table contains duplicate record_uid rows.")
  }

  tcols <- grep("^t[0-9]+$", names(x), value = TRUE)
  for (nm in tcols) {
    x[, (nm) := suppressWarnings(as.numeric(get(nm)))]
  }

  if ("expected_npeaks" %in% names(x)) {
    x[, expected_npeaks := suppressWarnings(as.integer(expected_npeaks))]
  }

  if ("exclude" %in% names(x)) {
    ex <- tolower(trimws(as.character(x$exclude)))
    x[, exclude := fifelse(
      is.na(ex) | ex == "",
      NA,
      ex %in% c("true", "t", "1", "yes", "y")
    )]
  }

  if ("note" %in% names(x)) {
    x[, note := as.character(note)]
  }

  x[]
}

lookup_rows_with_info <- function(lookup) {
  if (is.null(lookup) || !nrow(lookup)) return(data.table())

  x   <- as.data.table(copy(lookup))
  has <- rep(FALSE, nrow(x))

  tcols <- grep("^t[0-9]+$", names(x), value = TRUE)
  for (nm in tcols) {
    z   <- suppressWarnings(as.numeric(x[[nm]]))
    has <- has | is.finite(z)
  }

  if ("expected_npeaks" %in% names(x)) {
    z   <- suppressWarnings(as.numeric(x$expected_npeaks))
    has <- has | is.finite(z)
  }

  if ("exclude" %in% names(x)) {
    ## TRUE is an actionable exclusion. FALSE is treated as the blank/default
    ## state for backward compatibility with older templates that filled the
    ## entire column with FALSE. update_with_lookup() separately recognizes an
    ## explicit FALSE when it is used to un-exclude an already excluded run.
    has <- has | (!is.na(x$exclude) & x$exclude)
  }

  if ("note" %in% names(x)) {
    z   <- trimws(as.character(x$note))
    has <- has | (!is.na(z) & nzchar(z))
  }

  x[has]
}

.lookup_values <- function(lr) {
  out <- list(
    excluded        = FALSE,
    exclude_set     = FALSE,
    note            = "",
    expected_npeaks = NA_integer_,
    t0              = NA_real_,
    analyte_hints   = numeric()
  )

  if (is.null(lr) || !nrow(lr)) return(out)

  if ("exclude" %in% names(lr) && !is.na(lr$exclude[1L])) {
    out$exclude_set <- TRUE
    out$excluded    <- isTRUE(as.logical(lr$exclude[1L]))
  }

  if ("note" %in% names(lr) && !is.na(lr$note[1L])) {
    out$note <- trimws(as.character(lr$note[1L]))
  }

  if ("expected_npeaks" %in% names(lr) && is.finite(lr$expected_npeaks[1L])) {
    out$expected_npeaks <- as.integer(lr$expected_npeaks[1L])
  }

  if ("t0" %in% names(lr) && is.finite(lr$t0[1L])) {
    out$t0 <- as.numeric(lr$t0[1L])
  }

  tcols <- grep("^t[1-9][0-9]*$", names(lr), value = TRUE)
  if (length(tcols)) {
    tcols <- tcols[order(as.integer(sub("^t", "", tcols)))]
    z     <- suppressWarnings(as.numeric(unlist(lr[, ..tcols], use.names = FALSE)))
    out$analyte_hints <- z[is.finite(z)]
  }

  out
}

# ---------- Baseline preprocessing -------------------------------------------

auto_mile_mollify_hw <- function(
    hw,
    fraction          = 0.10,
    max_hw            = 15L,
    min_hw_to_mollify = 8L
  ) {

  hw <- suppressWarnings(as.integer(round(hw)))
  if (!is.finite(hw) || hw < min_hw_to_mollify) return(0L)

  max_hw <- suppressWarnings(as.integer(max_hw))
  if (!is.finite(max_hw) || max_hw < 1L) return(0L)

  as.integer(max(1L, min(max_hw, round(fraction * hw))))
}

preprocess_chromatogram <- function(
    d,
    baseline_method          = c("flatfit", "mile"),
    flatfit_args             = list(),
    mile_args                = list(),
    mile_auto_mollify        = TRUE,
    mile_mollify_fraction    = 0.10,
    mile_mollify_max_hw      = 15L,
    mile_mollify_min_base_hw = 8L
  ) {

  baseline_method <- match.arg(baseline_method)

  d <- copy(d)[order(time_min)]
  
  if (baseline_method == "flatfit") {
    
    args   <- flatfit_args
    args$y <- d$intensity
    args$x <- d$time_min
    
    mf <- do.call(flatfit, args)
    
    if (length(mf$baseline) != nrow(d) || length(mf$corrected) != nrow(d)) {
      stop("mile() returned vectors of unexpected length for: ", d$filename[1L])
    }
    
    d[, `:=`(
      baseline            = as.numeric(mf$baseline),
      intensity_corrected = as.numeric(mf$corrected),
      flatfit_weights     = as.numeric(mf$weights),
      flatfit_lambda      = as.numeric(mf$lambda),
      flatfit_smoothness  = as.numeric(mf$smoothness),
      flatfit_sg_hw       = as.numeric(mf$sg_hw),
      flatfit_iterations  = as.numeric(mf$iterations)
    )]
    
    attr(d, "flatfit_result") <- mf
    
  } else {
    
    args   <- mile_args
    args$y <- d$intensity
    args$x <- d$time_min
    if (is.null(args$direction)) args$direction <- "up"
    
    ## A small final mollifier is useful for removing knot-scale ripple from the
    ## estimated baseline, but it should remain much narrower than a genuine peak.
    ## An explicit mile_args$mollify_hw always overrides this rule.
    if (is.null(args$mollify_hw) && isTRUE(mile_auto_mollify)) {
      if (is.null(args$hw)) {
        z_for_hw <- if (identical(args$direction, "down")) -d$intensity else d$intensity
        args$hw <- estimate_hw(z_for_hw)
      }
      
      args$mollify_hw <- auto_mile_mollify_hw(
        hw                = args$hw,
        fraction          = mile_mollify_fraction,
        max_hw            = mile_mollify_max_hw,
        min_hw_to_mollify = mile_mollify_min_base_hw
      )
    }
    
    ## Preserve the original MILE default when automatic mollification is disabled
    ## and no explicit value was supplied.
    if (is.null(args$mollify_hw)) args$mollify_hw <- 0L
    
    mf <- do.call(mile, args)
    
    if (length(mf$baseline) != nrow(d) || length(mf$corrected) != nrow(d)) {
      stop("mile() returned vectors of unexpected length for: ", d$filename[1L])
    }
    
    d[, `:=`(
      baseline            = as.numeric(mf$baseline),
      intensity_corrected = as.numeric(mf$corrected),
      mile_sigma          = as.numeric(mf$sigma),
      mile_hw             = as.integer(mf$hw),
      mile_smooth_hw      = as.integer(mf$smooth_hw),
      mile_mollify_hw     = as.integer(args$mollify_hw),
      mile_iterations     = as.integer(mf$iterations)
    )]
    
    attr(d, "mile_result") <- mf
  }
  
  d[]
}

# ---------- Peak-selection helpers -------------------------------------------

interp_crossing <- function(x1, y1, x2, y2, target) {
  if (!all(is.finite(c(x1, y1, x2, y2, target)))) return(NA_real_)
  if (y2 == y1) return(mean(c(x1, x2)))
  x1 + (target - y1) * (x2 - x1) / (y2 - y1)
}

half_crossings <- function(t, y, idx, left_bound = 1L, right_bound = length(y)) {
  half <- y[idx] / 2

  left_x <- right_x <- NA_real_

  if (idx > left_bound) {
    jj <- which(y[left_bound:idx] <= half)
    if (length(jj)) {
      j <- left_bound + tail(jj, 1L) - 1L
      if (j < idx) {
        left_x <- interp_crossing(t[j], y[j], t[j + 1L], y[j + 1L], half)
      }
    }
  }

  if (idx < right_bound) {
    jj <- which(y[idx:right_bound] <= half)
    if (length(jj)) {
      j <- idx + jj[1L] - 1L
      if (j > idx) {
        right_x <- interp_crossing(t[j - 1L], y[j - 1L], t[j], y[j], half)
      }
    }
  }

  list(
    half  = half,
    left  = left_x,
    right = right_x,
    fwhm  = if (is.finite(left_x) && is.finite(right_x)) right_x - left_x else NA_real_
  )
}

trapz_base <- function(x, y) {
  if (length(x) < 2L) return(NA_real_)
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2, na.rm = TRUE)
}

snap_to_local_apex <- function(proc, hint_min, window_min = 0.10) {
  ii <- which(abs(proc$time_min - hint_min) <= window_min)
  if (!length(ii)) return(NA_integer_)
  ii[which.max(proc$intensity_corrected[ii])]
}

nearest_picker_peak <- function(pk, hint_min, window_min = 0.10) {
  if (is.null(pk) || !nrow(pk)) return(NA_integer_)
  ii <- which(abs(pk$rt - hint_min) <= window_min)
  if (!length(ii)) return(NA_integer_)

  dd   <- abs(pk$rt[ii] - hint_min)
  best <- ii[dd == min(dd)]
  if (length(best) > 1L) {
    best <- best[which.max(pk$height[best])]
  }
  as.integer(best[1L])
}

suppress_close_picker_peaks <- function(
    pk,
    min_peak_distance_min = 0.05,
    max_peaks             = Inf
  ) {

  if (!nrow(pk)) return(pk)

  ord    <- order(pk$snr, pk$height, decreasing = TRUE)
  chosen <- integer()

  for (j in ord) {
    if (!length(chosen) ||
        all(abs(pk$rt[j] - pk$rt[chosen]) >= min_peak_distance_min)) {
      chosen <- c(chosen, j)
      if (length(chosen) >= max_peaks) break
    }
  }

  pk[sort(chosen), ][order(rt)]
}

standardize_picker_peak <- function(
    proc,
    row,
    peak_type,
    peak_index,
    peak_source,
    noise_sd = NA_real_
  ) {

  use_deconv <- "deconv_accepted" %in% names(row) && isTRUE(row$deconv_accepted[1L])

  if (use_deconv) {
    rt  <- as.numeric(row$deconv_rt[1L])
    idx <- which.min(abs(proc$time_min - rt))
    baseline_rt <- as.numeric(stats::approx(
      proc$time_min,
      proc$baseline,
      xout = rt,
      rule = 2
    )$y)

    height <- as.numeric(row$deconv_height[1L])
    fwhm   <- as.numeric(row$deconv_fwhm[1L])

    return(data.table(
      idx                  = as.integer(idx),
      peak_type            = peak_type,
      peak_index           = as.integer(peak_index),
      peak_source          = peak_source,
      tR_min               = rt,
      apex_intensity       = baseline_rt + height,
      height               = height,
      prominence           = as.numeric(row$prominence[1L]),
      area                 = as.numeric(row$deconv_area[1L]),
      SNR                  = height / max(noise_sd, .Machine$double.eps),
      left_bound_min       = proc$time_min[as.integer(row$left[1L])],
      right_bound_min      = proc$time_min[as.integer(row$right[1L])],
      fwhm_left_min        = as.numeric(row$deconv_fwhm_left[1L]),
      fwhm_right_min       = as.numeric(row$deconv_fwhm_right[1L]),
      fwhm_min             = fwhm,
      fwhm_method          = "EGH_deconvolution",
      half_level_corrected = height / 2,
      half_level_raw       = baseline_rt + height / 2,
      fwhm_clean           = is.finite(fwhm) && isTRUE(row$deconv_converged[1L]),
      overlap_flag         = FALSE,
      edge_peak            = FALSE,
      saturated            = isTRUE(row$saturated[1L]),
      picker_type          = as.character(row$type[1L]),
      picker_cluster       = as.integer(row$cluster[1L]),
      picker_resolution    = as.numeric(row$resolution[1L]),
      clustered            = TRUE,
      deconvolved          = TRUE,
      deconv_sigma         = as.numeric(row$deconv_sigma[1L]),
      deconv_tau           = as.numeric(row$deconv_tau[1L]),
      deconv_k_detected    = as.integer(row$deconv_k_detected[1L]),
      deconv_k_fitted      = as.integer(row$deconv_k_fitted[1L])
    ))
  }

  idx <- as.integer(row$idx[1L])

  data.table(
    idx                  = idx,
    peak_type            = peak_type,
    peak_index           = as.integer(peak_index),
    peak_source          = peak_source,
    tR_min               = as.numeric(row$rt[1L]),
    apex_intensity       = proc$intensity[idx],
    height               = as.numeric(row$height[1L]),
    prominence           = as.numeric(row$prominence[1L]),
    area                 = as.numeric(row$area[1L]),
    SNR                  = as.numeric(row$snr[1L]),
    left_bound_min       = proc$time_min[as.integer(row$left[1L])],
    right_bound_min      = proc$time_min[as.integer(row$right[1L])],
    fwhm_left_min        = as.numeric(row$fwhm_left[1L]),
    fwhm_right_min       = as.numeric(row$fwhm_right[1L]),
    fwhm_min             = as.numeric(row$fwhm[1L]),
    fwhm_method          = "pick_peaks",
    half_level_corrected = as.numeric(row$height[1L]) / 2,
    half_level_raw       = proc$baseline[idx] + as.numeric(row$height[1L]) / 2,
    fwhm_clean           = is.finite(row$fwhm[1L]),
    overlap_flag         = identical(as.character(row$type[1L]), "shoulder"),
    edge_peak            = FALSE,
    saturated            = isTRUE(row$saturated[1L]),
    picker_type          = as.character(row$type[1L]),
    picker_cluster       = as.integer(row$cluster[1L]),
    picker_resolution    = as.numeric(row$resolution[1L]),
    clustered            = if ("deconv_attempted" %in% names(row)) isTRUE(row$deconv_attempted[1L]) else FALSE,
    deconvolved          = FALSE,
    deconv_sigma         = NA_real_,
    deconv_tau           = NA_real_,
    deconv_k_detected    = if ("deconv_k_detected" %in% names(row)) as.integer(row$deconv_k_detected[1L]) else NA_integer_,
    deconv_k_fitted      = if ("deconv_k_fitted" %in% names(row)) as.integer(row$deconv_k_fitted[1L]) else NA_integer_
  )
}

standardize_manual_peak <- function(
    proc,
    idx,
    selected_indices,
    peak_type,
    peak_index,
    peak_source,
    noise_sd
  ) {

  t   <- proc$time_min
  y   <- proc$intensity_corrected
  raw <- proc$intensity
  n   <- nrow(proc)

  selected_indices <- sort(unique(as.integer(selected_indices)))
  pos              <- match(idx, selected_indices)

  if (is.na(pos) || pos == 1L) {
    left_bound <- 1L
  } else {
    prev       <- selected_indices[pos - 1L]
    rr         <- prev:idx
    left_bound <- rr[which.min(y[rr])]
  }

  if (is.na(pos) || pos == length(selected_indices)) {
    right_bound <- n
  } else {
    nxt         <- selected_indices[pos + 1L]
    rr          <- idx:nxt
    right_bound <- rr[which.min(y[rr])]
  }

  hc <- half_crossings(t, y, idx, left_bound, right_bound)

  eps       <- max(abs(raw), na.rm = TRUE) * 1e-12 + .Machine$double.eps
  around    <- max(1L, idx - 2L):min(n, idx + 2L)
  saturated <- sum(abs(raw[around] - raw[idx]) <= eps) >= 3L

  left_overlap  <- !is.finite(hc$left)  && left_bound > 1L && y[left_bound] > hc$half
  right_overlap <- !is.finite(hc$right) && right_bound < n && y[right_bound] > hc$half

  data.table(
    idx                  = as.integer(idx),
    peak_type            = peak_type,
    peak_index           = as.integer(peak_index),
    peak_source          = peak_source,
    tR_min               = t[idx],
    apex_intensity       = raw[idx],
    height               = y[idx],
    prominence           = NA_real_,
    area                 = trapz_base(t[left_bound:right_bound], pmax(y[left_bound:right_bound], 0)),
    SNR                  = y[idx] / noise_sd,
    left_bound_min       = t[left_bound],
    right_bound_min      = t[right_bound],
    fwhm_left_min        = hc$left,
    fwhm_right_min       = hc$right,
    fwhm_min             = hc$fwhm,
    fwhm_method          = if (is.finite(hc$fwhm)) "manual_direct" else NA_character_,
    half_level_corrected = hc$half,
    half_level_raw       = proc$baseline[idx] + hc$half,
    fwhm_clean           = is.finite(hc$fwhm),
    overlap_flag         = left_overlap || right_overlap,
    edge_peak            = (!is.finite(hc$left) && left_bound == 1L) || (!is.finite(hc$right) && right_bound == n),
    saturated            = saturated,
    picker_type          = "manual",
    picker_cluster       = NA_integer_,
    picker_resolution    = NA_real_,
    clustered            = FALSE,
    deconvolved          = FALSE,
    deconv_sigma         = NA_real_,
    deconv_tau           = NA_real_,
    deconv_k_detected    = NA_integer_,
    deconv_k_fitted      = NA_integer_
  )
}

# ---------- Cluster deconvolution ---------------------------------------------

add_cluster_deconvolution <- function(
    proc,
    pk,
    analyte_selection,
    enabled            = TRUE,
    min_members        = 2L,
    pad                = 1.5,
    k_max              = 4L,
    min_snr            = 20,
    require_same_count = TRUE
  ) {

  pk <- copy(pk)

  deconv_cols <- list(
    deconv_attempted   = FALSE,
    deconv_accepted    = FALSE,
    deconv_rt          = NA_real_,
    deconv_height      = NA_real_,
    deconv_fwhm        = NA_real_,
    deconv_fwhm_left   = NA_real_,
    deconv_fwhm_right  = NA_real_,
    deconv_area        = NA_real_,
    deconv_sigma       = NA_real_,
    deconv_tau         = NA_real_,
    deconv_converged   = NA,
    deconv_k_detected  = NA_integer_,
    deconv_k_fitted    = NA_integer_,
    deconv_count_agree = NA
  )

  for (nm in names(deconv_cols)) {
    if (!nm %in% names(pk)) pk[, (nm) := deconv_cols[[nm]]]
  }

  info <- list(
    pk                    = pk,
    attempted_clusters    = integer(),
    accepted_clusters     = integer(),
    disagreement_clusters = integer(),
    failed_clusters       = integer(),
    fit                   = NULL,
    bic                   = NULL
  )

  if (!isTRUE(enabled) || !nrow(pk) || !nrow(analyte_selection)) return(info)

  picker_rows <- unique(analyte_selection[is.finite(picker_row), as.integer(picker_row)])
  picker_rows <- picker_rows[picker_rows >= 1L & picker_rows <= nrow(pk)]
  if (!length(picker_rows)) return(info)

  sel <- pk[picker_rows]
  clustered <- sel[, .N, by = cluster][
    is.finite(cluster) & N >= as.integer(min_members)
  ]
  if (!nrow(clustered)) return(info)

  cluster_ids <- clustered$cluster

  ## deconvolve_clusters() expects the original pick_peaks() schema.
  pk_fit <- as.data.frame(sel[cluster %in% cluster_ids])
  if (!"resolution" %in% names(pk_fit) && "picker_resolution" %in% names(pk_fit)) {
    pk_fit$resolution <- pk_fit$picker_resolution
  }

  fit <- try(
    deconvolve_clusters(
      x           = proc$time_min,
      y           = proc$intensity_corrected,
      pk          = pk_fit,
      min_members = as.integer(min_members),
      pad         = pad,
      k_max       = as.integer(k_max),
      min_snr     = min_snr
    ),
    silent = TRUE
  )

  info$attempted_clusters <- as.integer(cluster_ids)

  ## Mark every selected member of an attempted cluster even if fitting fails.
  for (cid in cluster_ids) {
    rows <- picker_rows[which(pk$cluster[picker_rows] == cid)]
    pk[rows, `:=`(
      deconv_attempted  = TRUE,
      deconv_k_detected = as.integer(length(rows))
    )]
  }

  if (inherits(fit, "try-error") || is.null(fit) || !nrow(fit)) {
    info$failed_clusters <- as.integer(cluster_ids)
    info$pk <- pk
    return(info)
  }

  fit_bic  <- attr(fit, "bic")
  fit      <- as.data.table(fit)
  info$fit <- fit
  info$bic <- fit_bic

  for (cid in cluster_ids) {
    rows <- picker_rows[pk$cluster[picker_rows] == cid]
    rows <- rows[order(pk$rt[rows])]
    ff   <- fit[cluster == cid][order(rt)]

    k_detected <- length(rows)
    k_fitted   <- nrow(ff)

    pk[rows, `:=`(
      deconv_k_detected  = as.integer(k_detected),
      deconv_k_fitted    = as.integer(k_fitted),
      deconv_count_agree = (k_detected == k_fitted)
    )]

    if (!k_fitted) {
      info$failed_clusters <- c(info$failed_clusters, as.integer(cid))
      next
    }

    if (isTRUE(require_same_count) && k_fitted != k_detected) {
      info$disagreement_clusters <- c(
        info$disagreement_clusters,
        as.integer(cid)
      )
      next
    }

    ## Automatic replacement is deliberately restricted to one-to-one detected
    ## versus fitted components. If different counts are ever allowed, that is a
    ## model-selection decision and should not be mapped silently by proximity.
    if (k_fitted != k_detected) {
      info$disagreement_clusters <- c(
        info$disagreement_clusters,
        as.integer(cid)
      )
      next
    }

    for (j in seq_along(rows)) {
      r <- rows[j]
      f <- ff[j]
      pk[r, `:=`(
        deconv_accepted   = TRUE,
        deconv_rt         = as.numeric(f$rt),
        deconv_height     = as.numeric(f$height),
        deconv_fwhm       = as.numeric(f$fwhm),
        deconv_fwhm_left  = as.numeric(f$fwhm_left),
        deconv_fwhm_right = as.numeric(f$fwhm_right),
        deconv_area       = as.numeric(f$area),
        deconv_sigma      = as.numeric(f$sigma),
        deconv_tau        = as.numeric(f$tau),
        deconv_converged  = isTRUE(f$converged)
      )]
    }

    info$accepted_clusters <- c(info$accepted_clusters, as.integer(cid))
  }

  info$attempted_clusters    <- unique(info$attempted_clusters)
  info$accepted_clusters     <- unique(info$accepted_clusters)
  info$disagreement_clusters <- unique(info$disagreement_clusters)
  info$failed_clusters       <- unique(info$failed_clusters)
  info$pk <- pk
  info
}

# ---------- One chromatogram --------------------------------------------------

analyze_one_chromatogram <- function(
    d,
    lookup                    = NULL,
    look_for_t0               = TRUE,
    t0_upper_limit            = NULL,
    t0_manual                 = NA_real_,
    analyte_lower_limit       = NULL,
    snr_min                   = 20,
    relative_height_min       = 0.25,
    fwer                      = 0.05,
    min_prom                  = NULL,
    nsim                      = 25L,
    do_shoulders              = FALSE,
    min_peak_distance_min     = 0.05,
    max_peaks                 = 6L,
    manual_window_min         = 0.10,
    baseline_method           = c("flatfit", "mile"),
    flatfit_args              = list(),
    mile_args                 = list(),
    mile_auto_mollify         = TRUE,
    mile_mollify_fraction     = 0.10,
    mile_mollify_max_hw       = 15L,
    mile_mollify_min_base_hw  = 8L,
    deconvolve_clustered      = TRUE,
    deconv_min_members        = 2L,
    deconv_pad                = 1.5,
    deconv_k_max              = 4L,
    deconv_min_snr            = NULL,
    deconv_require_same_count = TRUE
  ) {

  baseline_method <- match.arg(baseline_method)
  stopifnot(length(unique(d$record_uid)) == 1L)

  if (isTRUE(look_for_t0) &&
      (is.null(t0_upper_limit) || !is.finite(t0_upper_limit))) {
    stop("t0_upper_limit must be supplied when look_for_t0 = TRUE.")
  }

  if (length(t0_manual) != 1L) {
    stop("t0_manual must be a single numeric value or NA_real_.")
  }
  t0_manual <- suppressWarnings(as.numeric(t0_manual))
  if (is.finite(t0_manual) && t0_manual <= 0) {
    stop("t0_manual must be > 0 when supplied.")
  }

  if (!is.numeric(relative_height_min) || length(relative_height_min) != 1L ||
      !is.finite(relative_height_min) || relative_height_min < 0 ||
      relative_height_min > 1) {
    stop("relative_height_min must be a single value between 0 and 1.")
  }

  if (is.null(analyte_lower_limit)) {
    analyte_lower_limit <- if (isTRUE(look_for_t0)) t0_upper_limit else -Inf
  }

  proc <- preprocess_chromatogram(
    d,
    baseline_method          = baseline_method,
    flatfit_args             = flatfit_args,
    mile_args                = mile_args,
    mile_auto_mollify        = mile_auto_mollify,
    mile_mollify_fraction    = mile_mollify_fraction,
    mile_mollify_max_hw      = mile_mollify_max_hw,
    mile_mollify_min_base_hw = mile_mollify_min_base_hw
  )
  y    <- proc$intensity_corrected
  x    <- proc$time_min

  if (is.null(deconv_min_snr)) deconv_min_snr <- snr_min

  analyte_region <- which(x > analyte_lower_limit)

  ## Run the statistical detector with the absolute floor so a t0 marker is not
  ## lost merely because an analyte peak is extremely tall. The enantioseparation
  ## relative-height cutoff is applied only to analyte candidates afterwards.
  pk_raw <- pick_peaks(
    x            = x,
    y            = y,
    fwer         = fwer,
    min_prom     = min_prom,
    snr_min      = snr_min,
    nsim         = nsim,
    do_shoulders = do_shoulders
  )
  
  prom_cutoff <- attr(pk_raw, "prom_cutoff")
  noise_sd    <- as.numeric(attr(pk_raw, "noise_sd"))
  if (!is.finite(noise_sd) || noise_sd <= 0) {
    noise_sd <- .Machine$double.eps
  }
  pk <- as.data.table(pk_raw)

  if (nrow(pk)) {
    if (!"rt" %in% names(pk)) pk[, rt := x[idx]]
    if (!"snr" %in% names(pk)) pk[, snr := height / noise_sd]
  }

  ## Use the tallest statistically credible analyte candidate as the relative
  ## reference. If the picker found none, fall back to the tallest positive
  ## baseline-corrected sample in the analyte region. This avoids letting a
  ## single arbitrary raw point set the threshold whenever credible peaks exist.
  candidate_reference <- if (nrow(pk)) pk[rt > analyte_lower_limit] else data.table()
  max_analyte_signal <- if (nrow(candidate_reference)) {
    max(candidate_reference$height, na.rm = TRUE)
  } else if (length(analyte_region)) {
    max(pmax(y[analyte_region], 0), na.rm = TRUE)
  } else {
    NA_real_
  }
  
  analyte_height_floor <- max(
    snr_min * noise_sd,
    if (is.finite(max_analyte_signal)) relative_height_min * max_analyte_signal
    else -Inf
  )
  effective_analyte_snr_min <- analyte_height_floor / noise_sd
  relative_snr_min <- if (is.finite(max_analyte_signal)) {
    relative_height_min * max_analyte_signal / noise_sd
  } else {
    NA_real_
  }
  
  lv              <- .lookup_values(lookup)
  excluded        <- lv$excluded
  lookup_note     <- lv$note
  expected_npeaks <- lv$expected_npeaks
  lookup_t0       <- lv$t0
  manual_hints    <- lv$analyte_hints

  if (is.finite(lookup_t0) && lookup_t0 <= 0) {
    stop("Lookup t0 must be > 0 in file: ", proc$filename[1L])
  }

  manual_override <- is.finite(lookup_t0) ||
    length(manual_hints) > 0L ||
    is.finite(expected_npeaks) ||
    lv$exclude_set ||
    nzchar(lookup_note)

  ## ----- t0 marker / calculation value --------------------------------------

  t0_picker_row   <- NA_integer_
  t0_manual_idx   <- NA_integer_
  t0_source       <- "none"
  t0_value_min    <- NA_real_
  t0_marker_found <- FALSE

  if (isTRUE(look_for_t0)) {
    if (is.finite(lookup_t0)) {
      j <- nearest_picker_peak(pk, lookup_t0, manual_window_min)
      if (is.finite(j)) {
        t0_picker_row   <- j
        t0_marker_found <- TRUE
        t0_value_min    <- pk$rt[j]
        t0_source       <- "manual_hint_marker"
      } else {
        ii <- snap_to_local_apex(proc, lookup_t0, manual_window_min)
        if (is.finite(ii)) {
          t0_manual_idx   <- ii
          t0_marker_found <- TRUE
          t0_value_min    <- x[ii]
          t0_source       <- "manual_hint_marker"
        }
      }
    } else if (nrow(pk)) {
      c0 <- pk[rt <= t0_upper_limit]
      if (nrow(c0)) {
        j_local <- which.max(c0$height)
        t0_picker_row   <- which(pk$idx == c0$idx[j_local])[1L]
        t0_marker_found <- TRUE
        t0_value_min    <- c0$rt[j_local]
        t0_source       <- "automatic_marker"
      }
    }

    ## Optional exact fallback if the marker was expected but not found.
    if (!t0_marker_found && is.finite(t0_manual)) {
      t0_value_min <- t0_manual
      t0_source    <- "manual_value_fallback"
    }
  } else {
    ## With no marker, lookup t0 is an exact dead-time value rather than a peak
    ## hint. Per-file lookup takes precedence over a folder-level scalar.
    if (is.finite(lookup_t0)) {
      t0_value_min <- lookup_t0
      t0_source    <- "manual_value_lookup"
    } else if (is.finite(t0_manual)) {
      t0_value_min <- t0_manual
      t0_source    <- "manual_value"
    }
  }

  t0_available <- is.finite(t0_value_min) && t0_value_min > 0

  ## ----- analyte candidates --------------------------------------------------

  ca <- if (nrow(pk)) pk[rt > analyte_lower_limit] else data.table()

  if (nrow(ca)) {
    ca <- ca[snr >= effective_analyte_snr_min & height >= relative_height_min * max_analyte_signal]

    ca <- suppress_close_picker_peaks(
      ca,
      min_peak_distance_min = min_peak_distance_min,
      max_peaks             = max_peaks
    )
  }

  analyte_selection <- data.table(
    idx        = integer(),
    picker_row = integer(),
    source     = character()
  )

  if (length(manual_hints)) {
    for (hint in manual_hints) {
      j <- nearest_picker_peak(pk, hint, manual_window_min)

      if (is.finite(j) && pk$rt[j] > analyte_lower_limit) {
        analyte_selection <- rbind(
          analyte_selection,
          data.table(
            idx        = as.integer(pk$idx[j]),
            picker_row = as.integer(j),
            source     = "manual_hint"
          )
        )
      } else {
        ii <- snap_to_local_apex(proc, hint, manual_window_min)
        if (is.finite(ii) && x[ii] > analyte_lower_limit) {
          analyte_selection <- rbind(
            analyte_selection,
            data.table(
              idx        = as.integer(ii),
              picker_row = NA_integer_,
              source     = "manual_hint"
            )
          )
        }
      }
    }

    if (nrow(analyte_selection)) {
      analyte_selection <- unique(analyte_selection, by = "idx")
    }

    if (is.finite(expected_npeaks) &&
        nrow(analyte_selection) < expected_npeaks &&
        nrow(ca)) {

      ca_fill <- ca[order(-snr, -height)]
      for (j in seq_len(nrow(ca_fill))) {
        idx_j <- ca_fill$idx[j]
        rt_j  <- ca_fill$rt[j]

        if (!nrow(analyte_selection) ||
            all(abs(rt_j - x[analyte_selection$idx]) >= min_peak_distance_min)) {

          row_in_pk <- which(pk$idx == idx_j)[1L]
          analyte_selection <- rbind(
            analyte_selection,
            data.table(
              idx        = as.integer(idx_j),
              picker_row = as.integer(row_in_pk),
              source     = "automatic_fill"
            )
          )
        }

        if (nrow(analyte_selection) >= expected_npeaks) break
      }
    }
  } else if (nrow(ca)) {
    ca2 <- copy(ca)

    if (is.finite(expected_npeaks)) {
      ca2 <- ca2[order(-snr, -height)]
      ca2 <- ca2[seq_len(min(nrow(ca2), expected_npeaks))]
    }

    ca2 <- ca2[order(rt)]
    analyte_selection <- data.table(
      idx        = as.integer(ca2$idx),
      picker_row = as.integer(vapply(
        ca2$idx,
        function(ii) which(pk$idx == ii)[1L],
        integer(1)
      )),
      source     = "automatic"
    )
  }

  if (nrow(analyte_selection)) {
    analyte_selection <- unique(analyte_selection, by = "idx")
    setorder(analyte_selection, idx)

    if (nrow(analyte_selection) > max_peaks) {
      analyte_selection <- analyte_selection[seq_len(max_peaks)]
    }
  }

  ## ----- EGH refit for clustered, already-detected analyte peaks -----------

  deconv_info <- add_cluster_deconvolution(
    proc               = proc,
    pk                 = pk,
    analyte_selection  = analyte_selection,
    enabled            = deconvolve_clustered,
    min_members        = deconv_min_members,
    pad                = deconv_pad,
    k_max              = deconv_k_max,
    min_snr            = deconv_min_snr,
    require_same_count = deconv_require_same_count
  )
  pk <- deconv_info$pk

  ## ----- standardize selected peak rows ------------------------------------

  selected_indices <- analyte_selection$idx
  if (is.finite(t0_manual_idx)) selected_indices <- c(t0_manual_idx, selected_indices)
  if (is.finite(t0_picker_row)) selected_indices <- c(pk$idx[t0_picker_row], selected_indices)
  selected_indices <- sort(unique(as.integer(selected_indices)))

  peak_rows <- list()

  if (is.finite(t0_picker_row)) {
    peak_rows[[length(peak_rows) + 1L]] <- standardize_picker_peak(
      proc,
      pk[t0_picker_row],
      peak_type   = "t0",
      peak_index  = 0L,
      peak_source = t0_source,
      noise_sd    = noise_sd
    )
  } else if (is.finite(t0_manual_idx)) {
    peak_rows[[length(peak_rows) + 1L]] <- standardize_manual_peak(
      proc,
      idx              = t0_manual_idx,
      selected_indices = selected_indices,
      peak_type        = "t0",
      peak_index       = 0L,
      peak_source      = t0_source,
      noise_sd         = noise_sd
    )
  }

  if (nrow(analyte_selection)) {
    for (j in seq_len(nrow(analyte_selection))) {
      sel <- analyte_selection[j]

      if (is.finite(sel$picker_row)) {
        z <- standardize_picker_peak(
          proc,
          pk[sel$picker_row],
          peak_type   = "analyte",
          peak_index  = j,
          peak_source = sel$source,
          noise_sd    = noise_sd
        )
      } else {
        z <- standardize_manual_peak(
          proc,
          idx              = sel$idx,
          selected_indices = selected_indices,
          peak_type        = "analyte",
          peak_index       = j,
          peak_source      = sel$source,
          noise_sd         = noise_sd
        )
      }

      peak_rows[[length(peak_rows) + 1L]] <- z
    }
  }

  peaks <- if (length(peak_rows)) {
    rbindlist(peak_rows, use.names = TRUE, fill = TRUE)
  } else {
    data.table(
      idx                  = integer(),
      peak_type            = character(),
      peak_index           = integer(),
      peak_source          = character(),
      tR_min               = numeric(),
      apex_intensity       = numeric(),
      height               = numeric(),
      prominence           = numeric(),
      area                 = numeric(),
      SNR                  = numeric(),
      left_bound_min       = numeric(),
      right_bound_min      = numeric(),
      fwhm_left_min        = numeric(),
      fwhm_right_min       = numeric(),
      fwhm_min             = numeric(),
      fwhm_method          = character(),
      half_level_corrected = numeric(),
      half_level_raw       = numeric(),
      fwhm_clean           = logical(),
      overlap_flag         = logical(),
      edge_peak            = logical(),
      saturated            = logical(),
      picker_type          = character(),
      picker_cluster       = integer(),
      picker_resolution    = numeric(),
      clustered            = logical(),
      deconvolved          = logical(),
      deconv_sigma         = numeric(),
      deconv_tau           = numeric(),
      deconv_k_detected    = integer(),
      deconv_k_fitted      = integer(),
      k                    = numeric()
    )
  }

  if (nrow(peaks)) {
    peaks[, k := NA_real_]
    if (t0_available) {
      peaks[peak_type == "analyte", k := (tR_min - t0_value_min) / t0_value_min]
    }
  }

  meta_cols <- c(
    "record_uid",     "logical_key",  "date",        "folder_label",
    "source_folder",  "filename",     "file_hash",   "wavelength_nm",
    "experimentalID", "compoundname", "columnID",    "eluentID",
    "modifierID",     "temp_C",       "flow_mL_min", "metadata_raw"
  )

  if (nrow(peaks)) {
    meta <- unique(proc[, ..meta_cols])
    for (nm in meta_cols) peaks[, (nm) := meta[[nm]][1L]]
    setcolorder(peaks, c(meta_cols, setdiff(names(peaks), meta_cols)))
  }

  n_peaks <- nrow(peaks[peak_type == "analyte"])

  fwhm_review_any <- nrow(peaks[peak_type == "analyte"]) > 0L &&
    any(!is.finite(peaks[peak_type == "analyte", fwhm_min]))

  saturated_any <- nrow(peaks[peak_type == "analyte"]) > 0L &&
    any(peaks[peak_type == "analyte", saturated])

  clustered_any                        <- length(deconv_info$attempted_clusters) > 0L
  deconvolution_applied_any            <- length(deconv_info$accepted_clusters) > 0L
  deconvolution_count_disagreement_any <- length(deconv_info$disagreement_clusters) > 0L
  deconvolution_failed_any             <- length(deconv_info$failed_clusters) > 0L

  qc_status <- if (excluded) {
    "excluded"
  } else if (n_peaks == 0L) {
    "review_no_analyte_peak"
  } else if (isTRUE(look_for_t0) && !t0_marker_found && !t0_available) {
    "review_no_t0"
  } else if (deconvolution_count_disagreement_any) {
    "review_deconvolution_count"
  } else if (deconvolution_failed_any) {
    "review_deconvolution_fit"
  } else if (fwhm_review_any) {
    "review_fwhm"
  } else if (saturated_any) {
    "review_saturated"
  } else if (!isTRUE(look_for_t0) && !t0_available) {
    "ok_no_t0"
  } else if (manual_override || grepl("manual", t0_source, fixed = TRUE)) {
    "manual_guidance"
  } else {
    "ok"
  }

  run_meta <- unique(proc[, ..meta_cols])
  
  if (baseline_method == "flatfit") {
    baseline_info = list(
      flatfit_lambda     = proc$flatfit_lambda[1L],
      flatfit_smoothness = proc$flatfit_smoothness[1L],
      flatfit_sg_hw      = proc$flatfit_sg_hw[1L],
      flatfit_iterations = proc$flatfit_iterations[1L]
    )
  } else {
    baseline_info = list(
      mile_hw            = proc$mile_hw[1L],
      mile_smooth_hw     = proc$mile_smooth_hw[1L],
      mile_mollify_hw    = proc$mile_mollify_hw[1L],
      mile_iterations    = proc$mile_iterations[1L]
    )
  }

  run_summary <- run_meta[, .(
    record_uid,
    logical_key,
    date,
    folder_label,
    source_folder,
    filename,
    file_hash,
    wavelength_nm,
    experimentalID,
    compoundname,
    columnID,
    eluentID,
    modifierID,
    temp_C,
    flow_mL_min,
    metadata_raw,
    look_for_t0,
    t0_marker_found,
    t0_available,
    t0_value_min,
    t0_source,
    n_peaks,
    expected_npeaks,
    manual_override,
    excluded,
    fwhm_review_any,
    saturated_any,
    clustered_any,
    deconvolution_applied_any,
    deconvolution_count_disagreement_any,
    deconvolution_failed_any,
    peak_detection_ok = n_peaks >= 1L && !excluded && (!look_for_t0 || t0_marker_found || t0_available),
    noise_sd,
    max_analyte_signal,
    snr_min_absolute          = snr_min,
    relative_height_min,
    effective_analyte_snr_min,
    prominence_cutoff         = if (length(prom_cutoff)) as.numeric(prom_cutoff[1L]) else NA_real_,
    baseline_method,
    baseline_info             = list(baseline_info),
    qc_status,
    qc_note                   = lookup_note,
    peakfinder_version        = PEAKFINDER_VERSION,
    processed_at              = Sys.time()
  )]

  list(
    processed_trace = proc,
    candidates      = pk,
    peaks           = peaks,
    run             = run_summary,
    deconvolution   = deconv_info
  )
}

# ---------- QC plots ----------------------------------------------------------

plot_peak_qc <- function(result) {
  d   <- result$processed_trace
  pks <- result$peaks
  run <- result$run

  p <- ggplot(d, aes(time_min, intensity_corrected)) +
    geom_line(linewidth = 0.35) +
    geom_line(
      aes(y = baseline),
      linewidth = 0.3,
      linetype  = 2,
      color     = "deepskyblue"
    ) +
    labs(
      title    = d$filename[1L],
      subtitle = paste0(
        d$compoundname[1L], " | ", d$columnID[1L],           " | ",
        d$eluentID[1L],     " | ", d$modifierID[1L],         " | ",
        d$temp_C[1L],    " °C | ", d$flow_mL_min[1L], " mL/min | ",
        d$wavelength_nm[1L], " nm"
      ),
      x = "Time / min",
      y = "Intensity"
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

  if (nrow(pks)) {
    pp <- copy(pks)
    pp[, label := fifelse(
      peak_type == "t0",
      sprintf("t0 = %.3f", tR_min),
      sprintf(
        "tR%d = %.3f\nFWHM = %s%s",
        peak_index,
        tR_min,
        fifelse(is.finite(fwhm_min), sprintf("%.3f", fwhm_min), "NA"),
        fifelse(deconvolved, " [EGH]", "")
      )
    )]

    p <- p +
      geom_vline(
        data = pp,
        aes(xintercept = tR_min, color = peak_type),
        linewidth   = 0.35,
        linetype    = 3,
        inherit.aes = FALSE
      ) +
      geom_point(
        data = pp,
        aes(x = tR_min, y = apex_intensity, color = peak_type),
        size        = 1.0,
        inherit.aes = FALSE
      ) +
      geom_segment(
        data = pp[is.finite(fwhm_min)],
        aes(
          x     = fwhm_left_min,
          xend  = fwhm_right_min,
          y     = half_level_raw,
          yend  = half_level_raw,
          color = peak_type
        ),
        linewidth   = 1.0,
        inherit.aes = FALSE
      ) +
      geom_text(
        data = pp,
        aes(
          x     = tR_min,
          y     = apex_intensity,
          label = label,
          color = peak_type
        ),
        vjust         = -0.4,
        size          = 2.4,
        check_overlap = TRUE,
        inherit.aes   = FALSE
      ) +
      labs(color = "Peak")
  }

  ## If t0 is a supplied value rather than a visible marker, show it separately.
  if (isTRUE(run$t0_available[1L]) && !isTRUE(run$t0_marker_found[1L])) {
    p <- p +
      geom_vline(
        xintercept = run$t0_value_min[1L],
        linewidth  = 0.35,
        linetype   = 4
      ) +
      annotate(
        "text",
        x     = run$t0_value_min[1L],
        y     = Inf,
        label = sprintf("manual t0 = %.3f", run$t0_value_min[1L]),
        hjust = -0.05,
        vjust = 1.2,
        size  = 2.5
      )
  }

  p +
    annotate(
      "label",
      x     = Inf,
      y     = Inf,
      label = paste0(
        run$qc_status[1L],
        "\nSNR cutoff = ", sprintf("%.1f", run$effective_analyte_snr_min[1L])
      ),
      hjust = 1.05,
      vjust = 1.2,
      size  = 3
    )
}

save_peak_qc <- function(
    result,
    qc_dir,
    width  = 10,
    height = 4,
    dpi    = 160
  ) {

  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- tools::file_path_sans_ext(result$processed_trace$filename[1L])
  wl   <- result$processed_trace$wavelength_nm[1L]
  
  if (is.finite(wl)) {
    stem <- paste0(
      stem,
      "_",
      format(wl, trim = TRUE, scientific = FALSE),
      "nm"
    )
  }
  
  outfile <- file.path(qc_dir, paste0(stem, "_peakQC.png"))

  ggsave(
    outfile,
    plot   = plot_peak_qc(result),
    width  = width,
    height = height,
    dpi    = dpi
  )

  invisible(outfile)
}

# ---------- Wide ML/database table -------------------------------------------

snap_cross_wavelength_consensus <- function(peaks, rt_tolerance_min = 0.05) {
  
  x <- copy(peaks)
  
  if (!nrow(x)) {
    return(list(
      peaks               = x,
      peaks_for_retention = x,
      qc                  = data.table()
    ))
  }
  
  x[, `.row_id__` := .I]
  
  ## Preserve the original within-chromatogram assignment.
  x[, `:=`(
    peak_index_original                = peak_index,
    peak_index_consensus               = peak_index,
    cross_wavelength_consensus_member  = FALSE,
    cross_wavelength_consensus_support = NA_integer_,
    cross_wavelength_consensus_rt      = NA_real_,
    use_for_retention                  = TRUE
  )]
  
  groups <- x[
    peak_type == "analyte" &
      is.finite(tR_min) &
      !is.na(wavelength_nm),
    .(rows = list(.row_id__)),
    by = .(date, experimentalID)
  ]
  
  qc_list <- vector("list", nrow(groups))
  
  for (g in seq_len(nrow(groups))) {
    
    ii <- groups$rows[[g]]
    z  <- copy(x[.row_id__ %in% ii])
    
    n_wavelengths <- uniqueN(z$wavelength_nm)
    
    ## Nothing to align with only one wavelength.
    if (n_wavelengths < 2L) {
      qc_list[[g]] <- data.table(
        date                               = groups$date[g],
        experimentalID                     = groups$experimentalID[g],
        consensus_n_wavelengths            = n_wavelengths,
        consensus_required_support         = NA_integer_,
        consensus_n_peaks                  = NA_integer_,
        consensus_n_reindexed              = 0L,
        consensus_n_unmatched              = 0L,
        cross_wavelength_consensus_applied = FALSE,
        cross_wavelength_consensus_status  = "not_applicable"
      )
      next
    }
    
    majority_required <- floor(n_wavelengths / 2) + 1L
    
    ## Complete linkage is intentional:
    ## cutting at h = tolerance guarantees that the full spread of a cluster
    ## does not exceed rt_tolerance_min.
    if (nrow(z) == 1L) {
      cl <- 1L
    } else {
      cl <- stats::cutree(
        stats::hclust(
          stats::dist(z$tR_min),
          method = "complete"
        ),
        h = rt_tolerance_min
      )
    }
    
    z[, consensus_cluster := as.integer(cl)]
    
    cs <- z[
      ,
      .(
        consensus_rt         = stats::median(tR_min),
        rt_span              = max(tR_min) - min(tR_min),
        support              = uniqueN(wavelength_nm),
        duplicate_wavelength = anyDuplicated(wavelength_nm) > 0L
      ),
      by = consensus_cluster
    ]
    
    ## If two peaks from the same wavelength fall into the same tolerance
    ## cluster, their identity is intrinsically ambiguous at this tolerance.
    if (any(cs$duplicate_wavelength)) {
      
      qc_list[[g]] <- data.table(
        date                               = groups$date[g],
        experimentalID                     = groups$experimentalID[g],
        consensus_n_wavelengths            = n_wavelengths,
        consensus_required_support         = majority_required,
        consensus_n_peaks                  = 0L,
        consensus_n_reindexed              = 0L,
        consensus_n_unmatched              = 0L,
        cross_wavelength_consensus_applied = FALSE,
        cross_wavelength_consensus_status  = "ambiguous_no_change"
      )
      next
    }
    
    ## Only strict-majority clusters define trusted peak identities.
    maj <- cs[
      support >= majority_required &
        rt_span <= rt_tolerance_min
    ][
      order(consensus_rt)
    ]
    
    ## No majority -> preserve the original peak numbering.
    if (!nrow(maj)) {
      
      qc_list[[g]] <- data.table(
        date                               = groups$date[g],
        experimentalID                     = groups$experimentalID[g],
        consensus_n_wavelengths            = n_wavelengths,
        consensus_required_support         = majority_required,
        consensus_n_peaks                  = 0L,
        consensus_n_reindexed              = 0L,
        consensus_n_unmatched              = 0L,
        cross_wavelength_consensus_applied = FALSE,
        cross_wavelength_consensus_status  = "no_majority_no_change"
      )
      next
    }
    
    ## Consensus peak numbers are ordered by retention time.
    maj[, consensus_peak_index := seq_len(.N)]
    
    ## Once a reliable consensus exists, only majority-supported peaks are
    ## used to construct the aligned retention table.
    x[.row_id__ %in% ii, use_for_retention := FALSE]
    
    matched_ids <- integer()
    
    for (j in seq_len(nrow(maj))) {
      
      cid <- maj$consensus_cluster[j]
      ids <- z[consensus_cluster == cid, .row_id__]
      matched_ids <- c(matched_ids, ids)
      
      x[
        .row_id__ %in% ids,
        `:=`(
          peak_index                         = maj$consensus_peak_index[j],
          peak_index_consensus               = maj$consensus_peak_index[j],
          cross_wavelength_consensus_member  = TRUE,
          cross_wavelength_consensus_support = maj$support[j],
          cross_wavelength_consensus_rt      = maj$consensus_rt[j],
          use_for_retention                  = TRUE
        )
      ]
    }
    
    n_reindexed <- x[.row_id__ %in% matched_ids, sum(peak_index != peak_index_original)]
    n_unmatched <- x[.row_id__ %in% ii, sum(!use_for_retention)]
    
    qc_list[[g]] <- data.table(
      date                               = groups$date[g],
      experimentalID                     = groups$experimentalID[g],
      consensus_n_wavelengths            = n_wavelengths,
      consensus_required_support         = majority_required,
      consensus_n_peaks                  = nrow(maj),
      consensus_n_reindexed              = n_reindexed,
      consensus_n_unmatched              = n_unmatched,
      cross_wavelength_consensus_applied = TRUE,
      cross_wavelength_consensus_status  = if (n_reindexed > 0L) "reindexed" else "already_aligned"
    )
  }
  
  qc <- rbindlist(
    qc_list,
    use.names = TRUE,
    fill      = TRUE
  )
  
  x[, `.row_id__` := NULL]
  
  list(
    ## Complete audit copy: includes original and consensus indices.
    peaks = x[],
    
    ## This is what should feed make_retention_table().
    peaks_for_retention = x[use_for_retention == TRUE],
    
    qc = qc
  )
}

make_retention_table <- function(runs, peaks) {
  out <- copy(runs)

  ## t0_value_min is the source of truth for calculations. It may come either
  ## from a measured marker or from a manually supplied value.
  out[, t0 := as.numeric(t0_value_min)]

  if (nrow(peaks)) {
    t0_geom <- peaks[
      peak_type == "t0",
      .(
        record_uid,
        fwhm0        = fwhm_min,
        t0_SNR       = SNR,
        fwhm0_method = fwhm_method
      )
    ]

    if (nrow(t0_geom)) {
      out <- merge(out, t0_geom, by = "record_uid", all.x = TRUE, sort = FALSE)
    }
  }

  ## Keep the schema stable even when the whole batch has no visible t0 marker.
  if (!"fwhm0" %in% names(out))        out[, fwhm0 := NA_real_]
  if (!"t0_SNR" %in% names(out))       out[, t0_SNR := NA_real_]
  if (!"fwhm0_method" %in% names(out)) out[, fwhm0_method := NA_character_]

  an <- if (nrow(peaks)) peaks[peak_type == "analyte"] else data.table()
  if (!nrow(an)) return(out[])

  max_peak <- max(an$peak_index, na.rm = TRUE)

  for (j in seq_len(max_peak)) {
    tmp <- an[
      peak_index == j,
      .(
        record_uid,
        tR          = tR_min,
        fwhm        = fwhm_min,
        fwhm_method = fwhm_method,
        height      = height,
        area        = area,
        SNR         = SNR,
        prominence  = prominence,
        saturated   = saturated,
        peak_source = peak_source
      )
    ]

    setnames(
      tmp,
      setdiff(names(tmp), "record_uid"),
      paste0(setdiff(names(tmp), "record_uid"), j)
    )

    out <- merge(out, tmp, by = "record_uid", all.x = TRUE, sort = FALSE)

    tR_name <- paste0("tR", j)
    k_name  <- paste0("k", j)

    out[, (k_name) := fifelse(
      is.finite(t0) & t0 > 0 & is.finite(get(tR_name)),
      (get(tR_name) - t0) / t0,
      NA_real_
    )]
  }

  if (max_peak >= 2L) {
    for (j in seq_len(max_peak - 1L)) {
      j2 <- j + 1L

      t1 <- paste0("tR",   j)
      t2 <- paste0("tR",   j2)
      f1 <- paste0("fwhm", j)
      f2 <- paste0("fwhm", j2)
      K1 <- paste0("k",    j)
      K2 <- paste0("k",    j2)

      rs_name    <- paste0("Rs",    j, j2)
      alpha_name <- paste0("alpha", j, j2)

      out[, (rs_name) := fifelse(
        is.finite(get(t1)) & is.finite(get(t2)) &
          is.finite(get(f1)) & is.finite(get(f2)) &
          (get(f1) + get(f2)) > 0,
        1.18 * (get(t2) - get(t1)) / (get(f1) + get(f2)),
        NA_real_
      )]

      out[, (alpha_name) := fifelse(
        is.finite(get(K1)) & is.finite(get(K2)) & get(K1) > 0,
        get(K2) / get(K1),
        NA_real_
      )]
    }
  }

  out[]
}

spread_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  max(x) - min(x)
}

add_cross_wavelength_qc <- function(
    retention_table,
    peaks,
    rt_tolerance_min = 0.05
  ) {

  x <- copy(retention_table)

  exp_qc <- x[
    ,
    .(
      n_wavelengths = uniqueN(wavelength_nm),
      n_peaks_min   = if (all(is.na(n_peaks))) NA_integer_ else min(n_peaks, na.rm = TRUE),
      n_peaks_max   = if (all(is.na(n_peaks))) NA_integer_ else max(n_peaks, na.rm = TRUE),
      t0_spread_min = spread_safe(t0)
    ),
    by = .(date, experimentalID)
  ]

  if (nrow(peaks)) {
    tr <- peaks[
      peak_type == "analyte",
      .(tR_spread_min = spread_safe(tR_min)),
      by = .(date, experimentalID, peak_index)
    ]

    tr2 <- tr[
      ,
      .(
        max_tR_spread_min = {
          z <- tR_spread_min[is.finite(tR_spread_min)]
          if (length(z)) max(z) else NA_real_
        }
      ),
      by = .(date, experimentalID)
    ]

    exp_qc <- merge(
      exp_qc,
      tr2,
      by    = c("date", "experimentalID"),
      all.x = TRUE
    )
  } else {
    exp_qc[, max_tR_spread_min := NA_real_]
  }

  exp_qc[, cross_wavelength_peakcount_ok := n_peaks_min == n_peaks_max]
  exp_qc[
    ,
    cross_wavelength_rt_ok := fifelse(
      n_wavelengths < 2L,
      NA,
      cross_wavelength_peakcount_ok &
        (is.na(t0_spread_min) | t0_spread_min <= rt_tolerance_min) &
        (is.na(max_tR_spread_min) | max_tR_spread_min <= rt_tolerance_min)
    )
  ]

  merge(
    x,
    exp_qc,
    by    = c("date", "experimentalID"),
    all.x = TRUE,
    sort  = FALSE
  )
}

# ---------- Result assembly ---------------------------------------------------

assemble_analysis_result <- function(
    raw,
    results,
    path,
    settings,
    lookup_filename,
    qc_subdir,
    save_qc,
    cross_wavelength_tolerance_min,
    snap_to_cross_wavelength_consensus = TRUE
  ) {

  runs <- rbindlist(
    lapply(results, `[[`, "run"),
    use.names = TRUE,
    fill      = TRUE
  )

  peak_list <- lapply(results, `[[`, "peaks")
  peak_list <- peak_list[vapply(peak_list, nrow, integer(1)) > 0L]
  peaks <- if (length(peak_list)) {
    rbindlist(peak_list, use.names = TRUE, fill = TRUE)
  } else {
    data.table()
  }

  if (isTRUE(snap_to_cross_wavelength_consensus) && nrow(peaks)) {
    consensus <- snap_cross_wavelength_consensus(peaks, rt_tolerance_min = cross_wavelength_tolerance_min)
    
    peaks_consensus     <- consensus$peaks
    peaks_for_retention <- consensus$peaks_for_retention
    consensus_qc        <- consensus$qc
  } else {
    peaks_consensus     <- NULL
    peaks_for_retention <- peaks
    consensus_qc        <- data.table()
  }
  
  retention <- make_retention_table(runs, peaks_for_retention)
  
  retention <- add_cross_wavelength_qc(
    retention,
    peaks_for_retention,
    rt_tolerance_min = cross_wavelength_tolerance_min
  )
  
  if (nrow(consensus_qc)) {
    
    retention <- merge(
      retention,
      consensus_qc,
      by    = c("date", "experimentalID"),
      all.x = TRUE,
      sort  = FALSE
    )
  }

  list(
    raw             = raw,
    peaks           = peaks,
    peaks_consensus = peaks_consensus,
    runs            = runs,
    retention       = retention,
    results         = results,
    path            = path,
    settings        = settings,
    lookup_filename = lookup_filename,
    qc_subdir       = qc_subdir,
    save_qc         = save_qc,
    
    cross_wavelength_tolerance_min     = cross_wavelength_tolerance_min,
    snap_to_cross_wavelength_consensus = snap_to_cross_wavelength_consensus
  )
}

# ---------- Parallel helper ---------------------------------------------------

.parallel_lapply <- function(
    X,
    FUN,
    n_workers = 1L,
    seed      = 1L,
    progress  = TRUE
  ) {

  n_workers <- as.integer(n_workers)
  if (!is.finite(n_workers) || n_workers < 1L) n_workers <- 1L

  n <- length(X)
  if (!n) return(vector("list", 0L))

  ## Sequential mode gets a dependency-free base-R progress bar.
  if (n_workers == 1L || n <= 1L) {
    if (!is.null(seed)) set.seed(seed)
    out <- vector("list", n)
    names(out) <- names(X)

    pb <- NULL
    if (isTRUE(progress) && interactive()) {
      pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
      on.exit(if (!is.null(pb)) close(pb), add = TRUE)
    }

    for (i in seq_len(n)) {
      out[[i]] <- FUN(X[[i]])
      if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    }
    return(out)
  }

  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    warning(
      "n_workers > 1 requested, but packages 'future' and 'future.apply' are ",
      "not installed. Falling back to sequential analysis."
    )
    return(.parallel_lapply(X, FUN, n_workers = 1L, seed = seed, progress = progress))
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  max_workers <- suppressWarnings(parallel::detectCores(logical = TRUE)) - 1L
  if (!is.finite(max_workers) || max_workers < 1L) max_workers <- 1L
  n_workers <- max(1L, min(n_workers, max_workers, n))

  future::plan(future::multisession, workers = n_workers)
  if (!is.null(seed)) set.seed(seed)

  run_future <- function(progressor = NULL) {
    future.apply::future_lapply(
      X,
      function(z) {
        data.table::setDTthreads(1L)
        ans <- FUN(z)
        if (!is.null(progressor)) progressor()
        ans
      },
      future.seed     = TRUE,
      future.packages = c("data.table", "stringr")
    )
  }

  ## progressr relays progress conditions from multisession workers back to the
  ## main R process. Without it, future_lapply has no reliable completion-order
  ## progress signal, so analysis still runs but without a live bar.
  if (isTRUE(progress) && requireNamespace("progressr", quietly = TRUE)) {
    old_progress_opt <- getOption("progressr.enable")
    on.exit(options(progressr.enable = old_progress_opt), add = TRUE)
    options(progressr.enable = TRUE)

    out <- progressr::with_progress({
      p <- progressr::progressor(steps = n)
      run_future(progressor = p)
    })
  } else {
    if (isTRUE(progress) && !requireNamespace("progressr", quietly = TRUE)) {
      message(
        "Parallel analysis is running without a live progress bar. Install ",
        "'progressr' to enable one: install.packages('progressr')."
      )
    }
    out <- run_future()
  }

  names(out) <- names(X)
  out
}

# ---------- Folder-level analysis --------------------------------------------

analyze_chrom_folder <- function(
    path,
    look_for_t0                        = TRUE,
    t0_upper_limit                     = NULL,
    t0_manual                          = NA_real_,
    analyte_lower_limit                = NULL,
    snr_min                            = 20,
    relative_height_min                = 0.25,
    fwer                               = 0.05,
    min_prom                           = NULL,
    nsim                               = 25L,
    do_shoulders                       = FALSE,
    min_peak_distance_min              = 0.05,
    max_peaks                          = 6L,
    manual_window_min                  = 0.10,
    baseline_method                    = c("flatfit", "mile"),
    flatfit_args                       = list(),
    mile_args                          = list(),
    mile_auto_mollify                  = TRUE,
    mile_mollify_fraction              = 0.10,
    mile_mollify_max_hw                = 15L,
    mile_mollify_min_base_hw           = 8L,
    deconvolve_clustered               = TRUE,
    deconv_min_members                 = 2L,
    deconv_pad                         = 1.5,
    deconv_k_max                       = 4L,
    deconv_min_snr                     = NULL,
    deconv_require_same_count          = TRUE,
    lookup_filename                    = "peakfinder_lookup.csv",
    create_lookup_template             = TRUE,
    lookup_template_max_peaks          = 6L,
    overwrite_lookup_template          = FALSE,
    save_qc                            = TRUE,
    qc_subdir                          = "peakfinder_QC",
    cross_wavelength_tolerance_min     = 0.05,
    snap_to_cross_wavelength_consensus = TRUE,
    n_workers                          = 1L,
    seed                               = 1L,
    progress                           = TRUE,
    reader_args                        = list()
  ) {

  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  raw  <- do.call(read_chrom_folder, c(list(path), reader_args))
  
  record_map <- unique(raw[, .(record_uid, filename, wavelength_nm)])
  bad_uid    <- record_map[, .N, by = record_uid][N > 1L]
  
  if (nrow(bad_uid)) {
    stop("A record_uid maps to more than one chromatographic record.")
  }

  ## Important: the initial analysis never consumes the lookup table. The lookup
  ## is only applied by update_with_lookup(). create_lookup_template therefore
  ## means "write a blank operator template for a later correction pass".
  if (isTRUE(create_lookup_template)) {
    write_peak_lookup_template(
      raw,
      path      = path,
      filename  = lookup_filename,
      max_peaks = lookup_template_max_peaks,
      overwrite = overwrite_lookup_template
    )
  }

  settings <- list(
    look_for_t0               = look_for_t0,
    t0_upper_limit            = t0_upper_limit,
    t0_manual                 = t0_manual,
    analyte_lower_limit       = analyte_lower_limit,
    snr_min                   = snr_min,
    relative_height_min       = relative_height_min,
    fwer                      = fwer,
    min_prom                  = min_prom,
    nsim                      = nsim,
    do_shoulders              = do_shoulders,
    min_peak_distance_min     = min_peak_distance_min,
    max_peaks                 = max_peaks,
    manual_window_min         = manual_window_min,
    baseline_method           = baseline_method,
    flatfit_args              = flatfit_args,
    mile_args                 = mile_args,
    mile_auto_mollify         = mile_auto_mollify,
    mile_mollify_fraction     = mile_mollify_fraction,
    mile_mollify_max_hw       = mile_mollify_max_hw,
    mile_mollify_min_base_hw  = mile_mollify_min_base_hw,
    deconvolve_clustered      = deconvolve_clustered,
    deconv_min_members        = deconv_min_members,
    deconv_pad                = deconv_pad,
    deconv_k_max              = deconv_k_max,
    deconv_min_snr            = deconv_min_snr,
    deconv_require_same_count = deconv_require_same_count
  )

  record_ids <- unique(raw$record_uid)
  
  chrom <- setNames(
    lapply(
      record_ids,
      function(id) {
        raw[record_uid == id]
      }
    ),
    record_ids
  )

  worker <- function(dd) {
    do.call(
      analyze_one_chromatogram,
      c(list(d = dd, lookup = NULL), settings)
    )
  }

  results <- .parallel_lapply(
    chrom,
    worker,
    n_workers = n_workers,
    seed      = seed,
    progress  = progress
  )
  
  names(results) <- record_ids
  
  if (isTRUE(save_qc)) {
    qc_dir <- file.path(path, qc_subdir)
    for (i in seq_along(results)) {
      message("Saving QC plot ", i, " / ", length(results), ": ", names(results)[i])
      save_peak_qc(results[[i]], qc_dir = qc_dir)
    }
  }

  assemble_analysis_result(
    raw                                = raw,
    results                            = results,
    path                               = path,
    settings                           = settings,
    lookup_filename                    = lookup_filename,
    qc_subdir                          = qc_subdir,
    save_qc                            = save_qc,
    cross_wavelength_tolerance_min     = cross_wavelength_tolerance_min,
    snap_to_cross_wavelength_consensus = snap_to_cross_wavelength_consensus
  )
}

# ---------- Incremental lookup update ----------------------------------------

update_with_lookup <- function(
    result,
    lookup_filename = result$lookup_filename,
    save_qc         = result$save_qc,
    qc_subdir       = result$qc_subdir,
    n_workers       = 1L,
    seed            = 1L,
    progress        = TRUE
  ) {

  if (is.null(result$settings) || is.null(result$raw) || is.null(result$results)) {
    stop(
      "result does not contain the settings/raw/results fields required for an ",
      "incremental lookup update. Re-run analyze_chrom_folder() with this version."
    )
  }

  path   <- result$path
  lookup <- read_peak_lookup(path, lookup_filename)
  active <- lookup_rows_with_info(lookup)

  ## Backward-compatible un-exclusion: old templates commonly contain FALSE in
  ## every row. A FALSE becomes actionable only if that run is currently marked
  ## excluded in the result being updated.
  if ("exclude" %in% names(lookup) && "excluded" %in% names(result$runs)) {
    currently_excluded <- result$runs[excluded == TRUE, unique(record_uid)]
    unexclude <- lookup[!is.na(exclude) & exclude == FALSE & record_uid %in% currently_excluded]
    if (nrow(unexclude)) {
      active <- unique(rbindlist(list(active, unexclude), use.names = TRUE, fill = TRUE), by = "record_uid")
    }
  }

  if (!nrow(active)) {
    message("Lookup contains no rows with operator-specified information; nothing to update.")
    return(result)
  }

  known   <- names(result$results)
  unknown <- setdiff(active$record_uid, known)
  if (length(unknown)) {
    warning(
      "Ignoring lookup rows for filenames not present in the analyzed result: ",
      paste(unknown, collapse = ", ")
    )
  }

  target <- intersect(active$record_uid, known)
  if (!length(target)) {
    message("No lookup rows correspond to files in this result; nothing to update.")
    return(result)
  }

  message(
    "Updating ", length(target), " chromatogram(s) out of ", length(known),
    " using lookup guidance."
  )

  tasks <- setNames(
    lapply(target, function(id) list(
      data   = result$raw[record_uid == id],
      lookup = active[record_uid == id][1L]
    )),
    target
  )

  settings <- result$settings

  worker <- function(task) {
    do.call(
      analyze_one_chromatogram,
      c(list(d = task$data, lookup = task$lookup), settings)
    )
  }

  updated <- .parallel_lapply(
    tasks,
    worker,
    n_workers = n_workers,
    seed      = seed,
    progress  = progress
  )
  names(updated) <- target

  for (id in target) {
    result$results[[id]] <- updated[[id]]
  }

  if (isTRUE(save_qc)) {
    qc_dir <- file.path(path, qc_subdir)
    for (id in target) {
      save_peak_qc(result$results[[id]], qc_dir = qc_dir)
    }
  }

  out <- assemble_analysis_result(
    raw                                = result$raw,
    results                            = result$results,
    path                               = path,
    settings                           = result$settings,
    lookup_filename                    = lookup_filename,
    qc_subdir                          = qc_subdir,
    save_qc                            = save_qc,
    cross_wavelength_tolerance_min     = result$cross_wavelength_tolerance_min,
    snap_to_cross_wavelength_consensus = result$snap_to_cross_wavelength_consensus
  )

  out$updated_from_lookup <- target
  out
}

# ---------- Optional inspection pass -----------------------------------------

inspect_chrom_folder <- function(
    path,
    baseline_method = c("flatfit", "mile"),
    flatfit_args    = list(),
    mile_args       = list(),
    inspect_subdir  = "peakfinder_inspect",
    width           = 10,
    height          = 4,
    dpi             = 160,
    reader_args     = list()
  ) {

  path  <- normalizePath(path, winslash = "/", mustWork = TRUE)
  raw   <- do.call(read_chrom_folder, c(list(path), reader_args))
  
  record_ids  <- unique(raw$record_uid)
  inspect_dir <- file.path(path, inspect_subdir)
  dir.create(inspect_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_along(record_ids)) {
    id <- record_ids[i]
    dd <- raw[record_uid == id]
    
    message("Plotting chromatogram ", i, " / ", length(record_ids), ": ", dd$filename[1L])

    proc <- preprocess_chromatogram(
      d               = dd,
      baseline_method = baseline_method,
      flatfit_args    = flatfit_args,
      mile_args       = mile_args
    )

    p <- ggplot(proc, aes(time_min, intensity)) +
      geom_line(linewidth = 0.35) +
      geom_line(
        aes(y = baseline),
        linewidth = 0.3,
        linetype  = 2,
        color     = "deepskyblue"
      ) +
      labs(
        title    = proc$filename[1L],
        subtitle = paste0(
          proc$compoundname[1L], " | ", proc$columnID[1L],           " | ",
          proc$eluentID[1L],     " | ", proc$modifierID[1L],         " | ",
          proc$temp_C[1L],    " °C | ", proc$flow_mL_min[1L], " mL/min | ",
          proc$wavelength_nm[1L], " nm"
        ),
        x = "Time / min",
        y = "Intensity"
      ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())

    stem <- tools::file_path_sans_ext(dd$filename[1L])
    wl   <- dd$wavelength_nm[1L]
    if (is.finite(wl)) {
      stem <- paste0(
        stem,
        "_",
        format(wl, trim = TRUE, scientific = FALSE),
        "nm"
      )
    }
    
    ggsave(
      file.path(inspect_dir, paste0(stem, "_plot.png")),
      plot   = p,
      width  = width,
      height = height,
      dpi    = dpi
    )
  }

  invisible(inspect_dir)
}

# ---------- Duplicate-aware database append ----------------------------------

append_retention_database <- function(
    new_rows,
    database_file = "retention_database.rds",
    duplicate_key = c("date", "experimentalID", "wavelength_nm")
  ) {

  new_rows <- as.data.table(copy(new_rows))

  missing_key <- setdiff(duplicate_key, names(new_rows))
  if (length(missing_key)) {
    stop("Missing duplicate-key columns: ", paste(missing_key, collapse = ", "))
  }

  if (file.exists(database_file)) {
    old <- as.data.table(readRDS(database_file))
    db  <- rbindlist(list(old, new_rows), use.names = TRUE, fill = TRUE)
  } else {
    db <- new_rows
  }

  db[, duplicate := .N > 1L, by = duplicate_key]

  if (all(c("file_hash", "wavelength_nm") %in% names(db))) {
    db[, exact_file_duplicate := .N > 1L, by = .(file_hash, wavelength_nm)]
  } else if ("file_hash" %in% names(db)) {
    db[, exact_file_duplicate := .N > 1L, by = file_hash]
  } else {
    db[, exact_file_duplicate := NA]
  }

  dir.create(dirname(database_file), recursive = TRUE, showWarnings = FALSE)

  tmp <- tempfile(
    pattern = ".retention_database_",
    tmpdir  = dirname(database_file),
    fileext = ".rds"
  )
  saveRDS(db, tmp)

  ## Verify readability before replacing the current database.
  invisible(readRDS(tmp))

  if (file.exists(database_file)) {
    backup <- paste0(database_file, ".bak")
    if (file.exists(backup)) unlink(backup)

    if (!file.rename(database_file, backup)) {
      unlink(tmp)
      stop("Could not create database backup before replacement.")
    }

    if (!file.rename(tmp, database_file)) {
      file.rename(backup, database_file)
      stop("Could not replace database; original database was restored.")
    }

    unlink(backup)
  } else {
    if (!file.rename(tmp, database_file)) {
      unlink(tmp)
      stop("Could not create database file.")
    }
  }

  db[]
}

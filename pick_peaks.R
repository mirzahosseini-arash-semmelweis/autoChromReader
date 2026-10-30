## ---------------------------------------------------------------------------
##  Robust peak picking for (baseline-corrected) chromatograms
##
##  [MODIFIED v0.6.0] Optional DENOISING stage (section "Denoising"), the
##  bi-directional EMG peak model replacing EGH (section "BEMG peak shape"),
##  a faster ridge linker, a curvature-valley test for CWT shoulders, and a
##  full record of every resolved hyperparameter (attr "hyperparameters").
##
##  Three independent detectors, deliberately, because no single statistic
##  finds every kind of peak -- preceded by an optional denoiser:
##
##  (0) OPTIONAL DENOISING (default "none"). Savitzky-Golay, Gaussian-weighted
##      moving average, moving average, LOESS, Whittaker smoother, smoothing
##      spline or Butterworth FFT low-pass. The single smoothing parameter is
##      tuned automatically per chromatogram by the residual-variance (RV) or
##      autocorrelation (ACF) criterion of Bosten, Van Broeck & Cabooter
##      (2023), or fixed by the user. The SAME filter, with the SAME parameter,
##      is applied to every noise surrogate in (2), so the false-positive
##      control stays exact after smoothing.
##
##  (1) TOPOGRAPHIC PROMINENCE, computed exactly via 1-D persistent homology.
##      Answers "how far do I have to descend from this apex before I can
##      climb to a higher one". Scale free: a 5000-count peak and a 50-count
##      peak are treated identically.
##
##  (2) A NOISE-CALIBRATED CUTOFF for (1), obtained by block-bootstrapping the
##      high-frequency residual of the chromatogram itself and asking how large
##      a prominence pure noise produces. This is the operational answer to
##      "how would it know what is true": a peak is real if noise alone would
##      not have made it. Controls family-wise error over the whole run.
##
##  (3) A RICKER (Mexican-hat) CWT WITH RIDGE LINKING, i.e. the Du, Kibbe & Lin
##      (2006) scheme used by MassSpecWavelet and, in spirit, by xcms::centWave.
##      This is what finds MERGED peaks: a shoulder has almost no prominence
##      (it never descends to a valley) but it has a strong, scale-persistent
##      negative second derivative. The CWT with a Ricker wavelet is exactly a
##      smoothed second derivative, so a shoulder lights up a ridge.
##
##  Peaks found only by (3) are labelled "shoulder"; adjacent peaks are then
##  grouped into clusters and, optionally, fitted with a sum of
##  bi-directional exponentially modified Gaussian (BEMG) profiles so that the
##  apex of a genuinely merged component can be recovered rather than guessed.
##
##  What 1-D detection cannot do (documented limitation): two components that
##  co-elute so closely that the summed profile has neither a valley, nor a
##  separate negative-curvature lobe, nor a misfit a single BEMG cannot absorb
##  are indistinguishable from one asymmetric peak in a single channel. Only
##  a second data dimension (DAD spectra -> peak purity, EFA, MCR-ALS; or MS)
##  can reveal them. See pick_peaks() for the empirical resolution limits.
##
##  References
##    Du P., Kibbe W.A., Lin S.M. (2006) Bioinformatics 22:2059-2065.
##    Tautenhahn R. et al. (2008) BMC Bioinformatics 9:504.     (centWave)
##    Bosten E., Van Broeck P., Cabooter D. (2023) J. Chromatogr. A 1709:464360
##    Vivo-Truyols G., Schoenmakers P.J. (2006) Anal. Chem. 78:4598-4608.
##    Liitiainen E., Corona F., Lendasse A. (2010) Neurocomputing 73:1446.
##    Reed W.J. (2006) normal-Laplace distribution; Derrick J.T. et al. (2025)
##      J. Chromatogr. B 1264:124721 (PeakClimber); Wahab M.F., Handlovic T.T.
##      (2025) J. Sep. Sci. 48:e70139.                          (BEMG)
## ---------------------------------------------------------------------------


## --- Helper functions --------------------------------------------------------

## convolve y with a kernel of odd length, reflect-padded. Short kernels:
## direct O(n * length(k)); [MODIFIED] long kernels (> 101 taps): FFT,
## identical to rounding error and O(n log n).
.convolve <- function(y, k) {
  hw  <- (length(k) - 1L) %/% 2L
  if (hw < 1L) return(y * sum(k))
  n   <- length(y)
  hw  <- min(hw, n - 1L)
  if (length(k) > 2L * hw + 1L) {                 # kernel longer than trace
    mid <- (length(k) + 1L) %/% 2L
    k   <- k[(mid - hw):(mid + hw)]
  }
  if (length(k) > 101L) return(.convolve_fft(y, list(k))[, 1L])
  yp  <- .pad_reflect(y, hw)
  out <- numeric(n)
  for (j in seq_along(k)) out <- out + k[j] * yp[j:(j + n - 1L)]
  out
}

## Same operation for several odd-length kernels at once, sharing one FFT of
## the reflect-padded signal (used by the CWT, where all scales are needed).
.convolve_fft <- function(y, kernels) {
  n     <- length(y)
  hws   <- vapply(kernels, function(k) (length(k) - 1L) %/% 2L, integer(1))
  hwmax <- min(n - 1L, max(hws))
  yp    <- .pad_reflect(y, hwmax)
  ## zero-extend to a 2,3,5-smooth length: R's fft() is O(N^2) for lengths
  ## with large prime factors. The extension lies beyond the reflect padding
  ## and never reaches the returned samples.
  N     <- stats::nextn(length(yp))
  yp    <- c(yp, numeric(N - length(yp)))
  Fy    <- stats::fft(yp)
  out   <- matrix(0, n, length(kernels))
  for (j in seq_along(kernels)) {
    k  <- kernels[[j]]; hw <- hws[j]
    if (hw > hwmax) {
      mid <- hw + 1L; k <- k[(mid - hwmax):(mid + hwmax)]; hw <- hwmax
    }
    kk <- numeric(N)
    ## correlation kernel laid out circularly around lag 0
    kk[1L] <- k[hw + 1L]
    if (hw >= 1L) {
      kk[2:(hw + 1L)]         <- k[(hw + 2L):(2L * hw + 1L)]
      kk[(N - hw + 1L):N]     <- k[1:hw]
    }
    r <- Re(stats::fft(Fy * Conj(stats::fft(kk)), inverse = TRUE)) / N
    out[, j] <- r[(hwmax + 1L):(hwmax + n)]
  }
  out
}
## [MODIFIED] Savitzky-Golay coefficients: polynomial order `ord`, derivative
## `d`. Built on the scaled abscissa u = z / hw: exact, and numerically stable
## for long windows (raw powers make crossprod(X) singular beyond ~150 points).
.sg_coef <- function(hw, ord = 3L, d = 0L) {
  u <- (-hw:hw) / hw
  X <- outer(u, 0:ord, "^")
  C <- solve(crossprod(X), t(X))
  factorial(d) * C[d + 1L, ] / hw^d
}

sg_filter <- function(y, hw, ord = 3L, d = 0L, dx = 1) {
  if (hw < ord) hw <- ord + 1L
  .convolve(y, .sg_coef(hw, ord, d)) / dx^d
}

#' Noise measured the way a chromatographer measures it: in peak-free windows.
#'
#' The trace is cut into consecutive, non-overlapping windows of `win` points
#' (default 201, capped at n/8 so that short traces still give >= 8 windows).
#' A straight line is removed from each window -- this cancels residual drift,
#' so the per-window sd measures noise only -- and the per-window sd is
#' computed. Windows that contain any part of a peak have a much larger sd, so
#' the noise level is taken as the LOWER QUARTILE (q = 0.25) of the per-window
#' sds. q = 0.25 is a fixed, data-independent choice: it is robust as long as
#' at least ~25% of the windows are peak-free (true for practically every
#' HPLC run; the median would break down once peaks cover half the run).
#' With 201-point windows the sd of each window is itself estimated to ~5%
#' relative precision. The pooled residual of the windows with sd <= 1.5 x
#' sigma is returned as the noise pool that the bootstrap in
#' calibrate_prominence() resamples (bootstrapping the whole residual would
#' smuggle peak curvature into the null and inflate every threshold).
noise_profile <- function(y, win = 201L, q = 0.25) {
  n     <- length(y)
  win   <- as.integer(min(win, max(21L, n %/% 20L)))
  st    <- seq(1L, n - win + 1L, by = win)
  t0    <- seq_len(win); tc <- t0 - mean(t0); den <- sum(tc^2)
  resid <- vector("list", length(st))
  sds   <- numeric(length(st))
  for (k in seq_along(st)) {
    seg        <- y[st[k]:(st[k] + win - 1L)]
    r          <- seg - (mean(seg) + tc * (sum(tc * seg) / den))
    resid[[k]] <- r
    sds[k]     <- stats::sd(r)
  }
  sigma <- as.numeric(stats::quantile(sds, q))
  ## keep only windows consistent with the quiet quantile -- using the
  ## median would pull peak-containing windows into the pool whenever peaks
  ## occupy much of the run, and every threshold downstream would inflate
  quiet <- which(sds <= 1.5 * sigma)
  if (length(quiet) < 3L) quiet <- order(sds)[seq_len(min(3L, length(sds)))]
  pool <- unlist(resid[quiet], use.names = FALSE) * sqrt(win / (win - 2))
  list(
    sd      = sigma,
    pool    = pool,
    win_sd  = sds,
    n_quiet = length(quiet),
    win     = win,
    q       = q
  )
}

## sd of a convolution coefficient under the empirical noise autocovariance --
## exact for any stationary noise, no simulation needed.
.coef_sd <- function(k, gam) {
  m <- length(k); L <- min(m - 1L, length(gam) - 1L)
  s <- gam[1] * sum(k^2)
  if (L >= 1L) for (lag in seq_len(L))
    s <- s + 2 * gam[lag + 1L] * sum(k[1:(m - lag)] * k[(1 + lag):m])
  if (is.finite(s) && s > 0) sqrt(s) else NA_real_
}

## Ricker / Mexican-hat wavelet, L2-normalised, at scale s (in points).
.ricker <- function(s) {
  hw <- max(1L, ceiling(4 * s))
  t  <- -hw:hw
  k  <- (2 / (sqrt(3 * s) * pi^0.25)) * (1 - (t / s)^2) * exp(-t^2 / (2 * s^2))
  k / sqrt(sum(k^2))
}


## --- Denoising [NEW in v0.6.0] -----------------------------------------------
##
## Every denoiser below has exactly ONE smoothing parameter, expressed in
## sampling points (or the method's native penalty), so a single automatic
## tuning rule serves them all:
##
##   method      parameter                       reference
##   "sg"        half-window h (cubic)           Savitzky & Golay 1964
##   "gwma"      Gaussian sd s (points)          Gaussian-weighted moving avg
##   "ma"        half-window h                   boxcar moving average
##   "loess"     half-window h (quadratic, tricube weights; Cleveland 1979)
##   "whittaker" penalty lambda (2nd differences; Eilers 2003)
##   "spline"    smoothing parameter `spar` of stats::smooth.spline (all knots)
##   "fft"       cut-off frequency fc (cycles/point), 4th-order Butterworth,
##               zero-phase (PeakClimber uses an FFT low-pass)
##
## Kernel methods (sg, gwma, ma, loess) are applied as exact convolution
## kernels on the uniform sampling grid (LOESS here is the non-robust local
## quadratic fit, identical to stats::loess(degree = 2, family = "gaussian",
## surface = "direct") in the interior).
##
## Automatic tuning (Bosten et al. 2023):
##   tune = "rv"  residual-variance method. The noise variance V is estimated
##                from the raw trace with the modified 1-nearest-neighbour
##                estimator median(|y_i - y_{i-1}| * |y_i - y_{i+1}|), which
##                ignores sharp peaks because it uses the median. That median
##                is 0.7087 sigma^2 for Gaussian white noise (computed for
##                this implementation, 2e7 draws), so V is divided by 0.7087
##                to make it consistent. The parameter on the grid whose mean
##                squared residual mean((y - f)^2) is closest to V is chosen:
##                "smooth until what you remove looks like the noise".
##   tune = "acf" autocorrelation method (Vivo-Truyols & Schoenmakers 2006 as
##                used by Bosten et al.): the lag-one autocorrelation statistic
##                1 - sum(diff(e)^2) / (2 sum(e^2)) of the residual e = y - f is
##                matched to the same statistic of the 10-90% trimmed first
##                difference of the raw trace.
##   tune = "none" use `param` exactly as given.
## Both rules are data-driven per chromatogram; nothing depends on the
## sampling rate or on the analyte, and the resolved parameter is recorded.
## Integer kernels (sg, ma, loess) are scanned on a geometric grid of
## half-widths 1..min(250, n/40) points; continuous parameters (gwma,
## whittaker, spline, fft) are located by bisection over the same range,
## because both criteria increase monotonically with the amount of smoothing.

DENOISE_METHODS <- c("none", "sg", "gwma", "ma", "loess", "whittaker",
                     "spline", "fft")

## Gaussian consistency constant of the Bosten median 1-NN estimator.
.RV_GAUSS_CONST <- 0.7087

.ma_kernel <- function(h) {
  h <- max(1L, as.integer(h)); rep(1, 2L * h + 1L) / (2L * h + 1L)
}

.gauss_kernel <- function(s) {
  s  <- max(0.3, s)
  hw <- max(1L, ceiling(4 * s))
  t  <- -hw:hw
  k  <- exp(-t^2 / (2 * s^2))
  k / sum(k)
}

.sg_kernel <- function(h, ord = 3L) {
  h <- max(as.integer(h), ord %/% 2L + 1L)
  .sg_coef(h, ord, 0L)
}

.loess_kernel <- function(h, degree = 2L) {
  h  <- max(as.integer(h), degree)
  j  <- -h:h
  w  <- (1 - abs(j / (h + 1))^3)^3
  X  <- outer(j / h, 0:degree, "^")
  Xw <- X * w
  as.numeric(solve(crossprod(Xw, X), t(Xw))[1L, ])
}

## zero-phase Butterworth low-pass of order `ord`, applied in the frequency
## domain to a reflect-padded copy (padding suppresses wrap-around).
.fft_lowpass <- function(y, fc, ord = 4L) {
  n   <- length(y)
  pad <- as.integer(min(n - 1L, max(16L, ceiling(3 / fc))))
  yp  <- .pad_reflect(y, pad)
  N   <- stats::nextn(length(yp))                 # fast FFT length
  yp  <- c(yp, numeric(N - length(yp)))
  f   <- c(0:floor(N / 2), -rev(seq_len(ceiling(N / 2) - 1L))) / N
  H   <- 1 / sqrt(1 + (abs(f) / fc)^(2 * ord))
  out <- Re(stats::fft(stats::fft(yp) * H, inverse = TRUE)) / N
  out[(pad + 1L):(pad + n)]
}

.spline_smooth <- function(y, spar) {
  stats::smooth.spline(seq_along(y), y, all.knots = TRUE, spar = spar,
                       keep.data = FALSE)$y
}

## Apply one denoiser with a fixed parameter. Plain switch() with direct calls
## (no function names as strings), so future's global detection ships every
## helper to parallel workers.
.dn_apply <- function(y, method, param, sg_order = 3L) {
  switch(method,
    none      = y,
    sg        = .convolve(y, .sg_kernel(param, sg_order)),
    ma        = .convolve(y, .ma_kernel(param)),
    gwma      = .convolve(y, .gauss_kernel(param)),
    loess     = .convolve(y, .loess_kernel(param)),
    whittaker = .ws_solve(y, rep(1, length(y)), param),
    spline    = .spline_smooth(y, param),
    fft       = .fft_lowpass(y, param),
    stop("Unknown denoise method: ", method)
  )
}

## A closure that applies the SAME filter to any other vector (e.g. noise
## surrogates). Built in its own small environment so it does not capture the
## trace it was tuned on.
.dn_make_applier <- function(method, param, sg_order = 3L) {
  force(method); force(param); force(sg_order)
  function(z) .dn_apply(z, method, param, sg_order)
}

.dn_grid <- function(method, n) {
  hmax <- max(4L, min(250L, n %/% 40L))
  switch(method,
    sg        = unique(as.integer(round(exp(seq(log(2),   log(hmax), length.out = 30L))))),
    ma        = unique(as.integer(round(exp(seq(log(1),   log(hmax), length.out = 30L))))),
    loess     = unique(as.integer(round(exp(seq(log(2),   log(hmax), length.out = 30L))))),
    gwma      = exp(seq(log(0.5), log(hmax / 2), length.out = 30L)),
    whittaker = 10^seq(-1, 9, by = 0.25),
    spline    = seq(-1.5, 1.5, by = 0.1),
    fft       = exp(seq(log(0.45), log(0.5 / hmax), length.out = 30L)),
    numeric(0)
  )
}

## Noise-variance target of the RV rule.
.rv_target <- function(y) {
  n <- length(y)
  if (n < 5L) return(stats::var(y))
  i <- 2:(n - 1L)
  stats::median(abs(y[i] - y[i - 1L]) * abs(y[i] - y[i + 1L])) / .RV_GAUSS_CONST
}

.acf_stat <- function(e) {
  s <- sum(e^2)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  1 - sum(diff(e)^2) / (2 * s)
}

.acf_target <- function(y) {
  e  <- diff(y)
  qq <- stats::quantile(e, c(0.10, 0.90), names = FALSE)
  e  <- e[e > qq[1L] & e < qq[2L]]
  e  <- e[e != 0]
  .acf_stat(e)
}

## FWHM (points) of a denoiser's impulse response: a sampling-independent
## description of "how much smoothing" was applied, for the provenance record.
.dn_kernel_fwhm <- function(method, param, sg_order = 3L, n_max = 4001L) {
  if (identical(method, "none")) return(0)
  L <- n_max
  e <- numeric(L); e[(L + 1L) %/% 2L] <- 1
  r <- try(.dn_apply(e, method, param, sg_order), silent = TRUE)
  if (inherits(r, "try-error") || !all(is.finite(r))) return(NA_real_)
  m <- max(r); if (m <= 0) return(NA_real_)
  i <- which.max(r)
  l <- i; while (l > 1L && r[l] > m / 2) l <- l - 1L
  u <- i; while (u < L  && r[u] > m / 2) u <- u + 1L
  as.numeric(u - l)
}

#' Denoise a (baseline-corrected) chromatogram.
#'
#' @param y        numeric signal (uniformly sampled).
#' @param method   one of DENOISE_METHODS.
#' @param param    smoothing parameter (see table above). NULL = tune.
#' @param tune     "acf" (default), "rv", or "none" (requires `param`).
#' @param grid     optional candidate parameter grid for tuning.
#' @param sg_order polynomial order for "sg" (2 and 3 give identical smoothing
#'                 coefficients; 3 is used for consistency with the rest).
#' @return list(y, method, param, tune, target, criterion, grid,
#'              kernel_fwhm_pts, removed_sd)
denoise_chromatogram <- function(
    y,
    method   = DENOISE_METHODS,
    param    = NULL,
    tune     = c("acf", "rv", "none"),
    grid     = NULL,
    sg_order = 3L
  ) {
  method <- match.arg(method)
  tune   <- match.arg(tune)
  n      <- length(y)

  if (method == "none" || n < 7L) {
    return(list(y = y, method = "none", param = NA_real_, tune = "none",
                target = NA_real_, criterion = NA_real_, grid = NULL,
                kernel_fwhm_pts = 0, removed_sd = 0))
  }

  if (!is.null(param)) tune <- "none"
  if (tune == "none" && is.null(param)) {
    stop("denoise_chromatogram(): tune = 'none' requires a numeric `param`.")
  }

  target <- NA_real_; crit <- NA_real_
  if (tune != "none") {
    crit_fun <- if (tune == "rv") {
      function(f) mean((y - f)^2)
    } else {
      function(f) .acf_stat(y - f)
    }
    target <- if (tune == "rv") .rv_target(y) else .acf_target(y)
    if (is.null(grid) && method %in% c("gwma", "whittaker", "spline", "fft")) {
      ## continuous parameter: both criteria grow monotonically with the
      ## amount of smoothing, so the match is found by bisection on a
      ## transformed scale s (larger s = more smoothing), 18 halvings of the
      ## same range the grid would cover.
      g0 <- .dn_grid(method, n)
      tr <- switch(method,
        gwma      = list(to = function(s) exp(s),  lo = log(min(g0)),  hi = log(max(g0))),
        whittaker = list(to = function(s) 10^s,    lo = log10(min(g0)), hi = log10(max(g0))),
        spline    = list(to = function(s) s,       lo = min(g0),        hi = max(g0)),
        fft       = list(to = function(s) exp(-s), lo = -log(max(g0)),  hi = -log(min(g0)))
      )
      ev <- function(s) {
        f <- .dn_apply(y, method, tr$to(s), sg_order)
        list(s = s, v = crit_fun(f) - target, f = f)
      }
      a <- ev(tr$lo); b <- ev(tr$hi)
      path <- list(a, b)
      if (is.finite(a$v) && a$v >= 0) {
        sel <- a
      } else if (is.finite(b$v) && b$v <= 0) {
        sel <- b
      } else {
        for (it in seq_len(18L)) {
          mid <- ev((a$s + b$s) / 2)
          path[[length(path) + 1L]] <- mid
          if (!is.finite(mid$v)) break
          if (mid$v < 0) a <- mid else b <- mid
        }
        sel <- if (abs(a$v) <= abs(b$v)) a else b
      }
      param <- tr$to(sel$s)
      f     <- sel$f
      grid  <- vapply(path, function(z) tr$to(z$s), numeric(1))
      crit  <- vapply(path, function(z) z$v + target, numeric(1))
    } else {
      if (is.null(grid)) grid <- .dn_grid(method, n)
      fits <- lapply(grid, function(p) .dn_apply(y, method, p, sg_order))
      crit <- vapply(fits, crit_fun, numeric(1))
      ok   <- is.finite(crit)
      if (!any(ok)) stop("denoise_chromatogram(): tuning criterion undefined.")
      best  <- which(ok)[which.min(abs(crit[ok] - target))]
      param <- grid[best]
      f     <- fits[[best]]
    }
  } else {
    f <- .dn_apply(y, method, param, sg_order)
  }

  list(
    y               = f,
    method          = method,
    param           = as.numeric(param),
    tune            = tune,
    target          = target,
    criterion       = crit,
    grid            = grid,
    kernel_fwhm_pts = .dn_kernel_fwhm(method, param, sg_order),
    removed_sd      = stats::sd(y - f),
    sg_order        = sg_order
  )
}


## --- 1-D persistent homology -------------------------------------------------

#' Topographic prominence of every local maximum, via union-find on the
#' superlevel-set filtration. Each local max is "born" at its own height and
#' "dies" when it merges into a taller component; prominence = birth - death.
#' Equivalent to scipy.signal.peak_prominences but computed globally in one
#' O(n log n) pass, and it also returns the merge partner, which tells you
#' which peak a shoulder belongs to.
#'
#' @return data.frame(idx, birth, death, prominence, left, right)
prominence <- function(y) {
  n_full <- length(y)
  if (n_full < 3L) {
    return(data.frame(idx = which.max(y), birth = max(y), death = min(y),
                      prominence = max(y) - min(y), left = 1L, right = n_full))
  }

  ## (a) keep endpoints and every point that is not strictly monotone
  d  <- diff(y)
  up <- d > 0; dn <- d < 0
  mono <- (up[-length(up)] & up[-1]) | (dn[-length(dn)] & dn[-1])
  keepv <- c(TRUE, !mono, TRUE)
  map <- which(keepv)
  z <- y[map]
  n <- length(z)

  ord <- order(z, decreasing = TRUE)

  parent  <- integer(n)
  root_pk <- integer(n)
  root_lo <- integer(n)
  root_hi <- integer(n)
  seen    <- logical(n)

  m   <- 0L
  cap <- n
  pk  <- integer(cap); bir <- numeric(cap); dea <- numeric(cap)
  lo  <- integer(cap); hi  <- integer(cap)

  find <- function(i) {
    r <- i
    while (parent[r] != r) r <- parent[r]
    while (parent[i] != r) { nx <- parent[i]; parent[i] <<- r; i <- nx }
    r
  }

  for (i in ord) {
    seen[i] <- TRUE
    l <- if (i > 1L && seen[i - 1L]) find(i - 1L) else 0L
    r <- if (i < n  && seen[i + 1L]) find(i + 1L) else 0L

    if (l == 0L && r == 0L) {
      parent[i] <- i; root_pk[i] <- i; root_lo[i] <- i; root_hi[i] <- i
    } else if (r == 0L) {
      parent[i] <- l; root_hi[l] <- i
    } else if (l == 0L) {
      parent[i] <- r; root_lo[r] <- i
    } else {
      if (z[root_pk[l]] >= z[root_pk[r]]) { old <- l; yng <- r } else { old <- r; yng <- l }
      m <- m + 1L
      pk[m]  <- root_pk[yng]
      bir[m] <- z[root_pk[yng]]
      dea[m] <- z[i]
      lo[m]  <- root_lo[yng]; hi[m] <- root_hi[yng]
      parent[yng] <- old; parent[i] <- old
      if (root_lo[yng] < root_lo[old]) root_lo[old] <- root_lo[yng]
      if (i          < root_lo[old])   root_lo[old] <- i
      if (root_hi[yng] > root_hi[old]) root_hi[old] <- root_hi[yng]
      if (i          > root_hi[old])   root_hi[old] <- i
    }
  }

  g <- find(ord[1])
  m <- m + 1L
  pk[m] <- root_pk[g]; bir[m] <- z[root_pk[g]]; dea[m] <- min(z)
  lo[m] <- 1L;         hi[m]  <- n

  ii <- seq_len(m)
  out <- data.frame(idx        = map[pk[ii]],
                    birth      = bir[ii],
                    death      = dea[ii],
                    prominence = bir[ii] - dea[ii],
                    left       = map[lo[ii]],
                    right      = map[hi[ii]])
  out[order(out$idx), ]
}

## --- Noise calibration -------------------------------------------------------

## Circular block bootstrap of the high-frequency residual, preserving the
## short-range autocorrelation that detector filtering puts into HPLC noise
## (ignore it and every threshold derived will be too low).
.surrogate <- function(r, block, n_out = length(r)) {
  n   <- length(r)
  nb  <- ceiling(n_out / block)
  st  <- sample.int(n, nb, replace = TRUE)
  idx <- as.vector(vapply(st, function(s) ((s + 0:(block - 1L) - 1L) %% n) + 1L,
                          integer(block)))
  r[idx[seq_len(n_out)]]
}

## Biased sample autocovariance at lags 0..L via FFT (same
## estimator as stats::acf(type = "covariance"), but O(n log n)).
.acov_fft <- function(s, L) {
  s <- s - mean(s); n <- length(s)
  N <- stats::nextn(2L * n)
  f <- stats::fft(c(s, numeric(N - n)))
  Re(stats::fft(Mod(f)^2, inverse = TRUE))[seq_len(L + 1L)] / (N * n)
}

## Autocorrelation length of the residual: first lag where acf drops below 1/e.
.acf_len <- function(r, max_lag = 50L) {
  a <- stats::acf(r, lag.max = max_lag, plot = FALSE)$acf[-1]
  w <- which(a < exp(-1))
  if (!length(w)) max_lag else w[1]
}

#' Largest prominence attributable to noise alone.
#'
#' Surrogate noise traces of full length are built by circular block bootstrap
#' of the quiet-window residual pool (block = 3 x the lag at which the
#' residual autocorrelation drops below 1/e, so detector smoothing is
#' preserved). The largest prominence in each surrogate is recorded; the
#' (1 - fwer) quantile of those maxima is the cutoff. Nothing here is a tuning
#' constant: the cutoff is in the units of the data and is re-derived for
#' every chromatogram from its own noise.
#'
#' @param fwer probability that a noise-only chromatogram yields ANY peak above
#'   the returned cutoff. 0.05 => at most a 1-in-20 chance of a single false
#'   positive anywhere in the run.
#' @param denoiser [NEW] optional function applied to every surrogate before
#'   its prominences are measured -- the same filter that was applied to the
#'   data, so the null distribution matches what the detector actually sees.
calibrate_prominence <- function(
    y,
    nsim     = 25L,
    fwer     = 0.05,
    np       = NULL,
    denoiser = NULL
  ) {
  if (is.null(np)) np <- noise_profile(y)
  r  <- np$pool
  b  <- max(2L, .acf_len(r) * 3L)
  n  <- length(y)
  lag_max <- min(200L, n - 1L)
  gsum    <- numeric(lag_max + 1L)
  mx <- vapply(seq_len(nsim), function(i) {
    s <- .surrogate(r, b, n)
    if (!is.null(denoiser)) s <- denoiser(s)
    gsum <<- gsum + .acov_fft(s, lag_max)
    max(prominence(s)$prominence)
  }, numeric(1))
  list(
    cutoff   = as.numeric(stats::quantile(mx, 1 - fwer)),
    sd       = np$sd,
    block    = b,
    acf_len  = b %/% 3L,
    null_max = mx,
    gamma    = gsum / nsim
  )
}

## --- CWT ridge detection -----------------------------------------------------

#' Ricker CWT plus ridge linking (Du et al. 2006).
#'
#' [MODIFIED] Same linking rule as before (coarse -> fine, nearest unused
#' local maximum within win_fac x scale, a ridge ends after gap_max missed
#' scales) but stored in flat vectors instead of growing lists, which makes it
#' ~10-50x faster on noise and lets the shoulder threshold use more surrogates.
#' Conflicts between ridges competing for one maximum go to the closer ridge.
#'
#' @param scales wavelet scales in POINTS. Should bracket sigma of the
#'   narrowest and widest expected peak.
#' @return list(W = n x nscale coefficient matrix, ridges = data.frame)
cwt_ridges <- function(y, scales, gap_max = 3L, win_fac = 2) {
  n  <- length(y)
  W  <- .convolve_fft(y, lapply(scales, .ricker))
  ns <- length(scales)

  cap   <- 4096L
  e_id  <- integer(cap); e_sc <- integer(cap); e_pos <- integer(cap); m <- 0L
  push  <- function(id, sc, pos) {
    k <- length(id)
    if (!k) return(invisible())
    if (m + k > length(e_id)) {
      grow  <- max(length(e_id), k) * 2L
      e_id  <<- c(e_id,  integer(grow))
      e_sc  <<- c(e_sc,  integer(grow))
      e_pos <<- c(e_pos, integer(grow))
    }
    ii <- (m + 1L):(m + k)
    e_id[ii] <<- id; e_sc[ii] <<- sc; e_pos[ii] <<- pos
    m <<- m + k
    invisible()
  }

  a_id <- integer(0); a_pos <- integer(0); a_gap <- integer(0); next_id <- 0L

  for (j in rev(seq_len(ns))) {
    w    <- W[, j]
    dw   <- diff(w)
    cand <- which(c(FALSE, dw > 0) & c(dw < 0, FALSE))
    nc   <- length(cand)
    used <- logical(nc)
    win  <- max(1L, ceiling(win_fac * scales[j]))
    linked <- logical(length(a_id))

    if (length(a_id) && nc) {
      k   <- findInterval(a_pos, cand)
      lo  <- pmax(k, 1L); hi <- pmin(k + 1L, nc)
      dlo <- abs(cand[lo] - a_pos); dhi <- abs(cand[hi] - a_pos)
      nn  <- ifelse(dhi < dlo, hi, lo)
      alt <- ifelse(dhi < dlo, lo, hi)
      dd  <- pmin(dlo, dhi)
      for (a in order(dd)) {
        if (dd[a] > win) break
        c1 <- nn[a]
        if (!used[c1]) {
          used[c1] <- TRUE; linked[a] <- TRUE; a_pos[a] <- cand[c1]
        } else {
          c2 <- alt[a]
          if (!used[c2] && abs(cand[c2] - a_pos[a]) <= win) {
            used[c2] <- TRUE; linked[a] <- TRUE; a_pos[a] <- cand[c2]
          }
        }
      }
    }
    if (length(a_id)) {
      push(a_id[linked], rep(j, sum(linked)), a_pos[linked])
      a_gap[linked]  <- 0L
      a_gap[!linked] <- a_gap[!linked] + 1L
      alive <- a_gap <= gap_max
      a_id <- a_id[alive]; a_pos <- a_pos[alive]; a_gap <- a_gap[alive]
    }
    if (any(!used)) {
      newc  <- cand[!used]
      ids   <- next_id + seq_along(newc); next_id <- next_id + length(newc)
      push(ids, rep(j, length(newc)), newc)
      a_id  <- c(a_id, ids); a_pos <- c(a_pos, newc); a_gap <- c(a_gap, integer(length(newc)))
    }
  }

  if (!m) return(list(W = W, ridges = NULL, scales = scales))
  e_id <- e_id[seq_len(m)]; e_sc <- e_sc[seq_len(m)]; e_pos <- e_pos[seq_len(m)]
  coef <- W[cbind(e_pos, e_sc)]

  o_best <- order(e_id, -coef)
  best   <- o_best[!duplicated(e_id[o_best])]
  o_fine <- order(e_id, e_sc)
  fine   <- o_fine[!duplicated(e_id[o_fine])]
  len    <- tabulate(e_id, nbins = next_id)

  rid <- data.frame(
    idx      = e_pos[best],
    idx_fine = e_pos[fine][match(e_id[best], e_id[fine])],
    scale    = scales[e_sc[best]],
    len      = len[e_id[best]],
    coef     = coef[best]
  )
  list(
    W      = W,
    ridges = rid[order(rid$idx), ],
    scales = scales
  )
}

## --- BEMG peak shape [NEW in v0.6.0; replaces EGH] ---------------------------
##
## Bi-directional exponentially modified Gaussian: a Gaussian (mu, sigma)
## convolved with an asymmetric two-sided exponential (Laplace) kernel with
## left time constant tau_left (fronting) and right time constant tau_right
## (tailing). Equivalently, a mixture of a fronting EMG and a tailing EMG that
## share mu and sigma, with weights tau_l / (tau_l + tau_r) and
## tau_r / (tau_l + tau_r). Unit-area form:
##
##   f(t) = 1 / (tau_l + tau_r) * [
##          exp( (mu - t)/tau_r + sigma^2/(2 tau_r^2) ) * Phi( (t - mu)/sigma - sigma/tau_r )
##        + exp( (t - mu)/tau_l + sigma^2/(2 tau_l^2) ) * Phi( (mu - t)/sigma - sigma/tau_l ) ]
##
## tau_l -> 0 gives the classic tailing EMG, tau_r -> 0 a fronting EMG, both
## -> 0 the Gaussian. It is log-concave, hence unimodal. Each exponential
## factor and log Phi are summed in LOG space (pnorm(..., log.p = TRUE)),
## which is the numerically stable evaluation: the naive product overflows
## (exp -> Inf times Phi -> 0) as soon as sigma/tau exceeds ~25.
## Because the model is written in unit-area form, the linear amplitude is
## the peak AREA, and apex height / apex time / FWHM are obtained numerically.

.bemg_unit <- function(t, mu, sg, tl, tr) {
  z  <- (t - mu) / sg
  lr <- (mu - t) / tr + sg^2 / (2 * tr^2) + stats::pnorm(z - sg / tr,  log.p = TRUE)
  ll <- (t - mu) / tl + sg^2 / (2 * tl^2) + stats::pnorm(-z - sg / tl, log.p = TRUE)
  (exp(lr) + exp(ll)) / (tl + tr)
}

#' Bi-directional EMG with total area `A`.
bemg <- function(t, A, mu, sigma, tau_left, tau_right) {
  A * .bemg_unit(t, mu, sigma, tau_left, tau_right)
}

#' Apex position, apex height (per unit area) and width at `frac` height of a
#' BEMG component. Root-bracketing is exact because the model is unimodal.
bemg_shape <- function(mu, sigma, tau_left, tau_right, frac = 0.5) {
  f   <- function(t) .bemg_unit(t, mu, sigma, tau_left, tau_right)
  sp  <- 4 * sigma
  op  <- stats::optimize(f, c(mu - sp - 2 * tau_left, mu + sp + 2 * tau_right),
                         maximum = TRUE, tol = 1e-10 * max(1, abs(mu)) + 1e-12)
  apx <- op$maximum; fm <- op$objective
  lvl <- frac * fm
  g   <- function(t) f(t) - lvl
  lo  <- apx - (6 * sigma + 12 * tau_left)
  hi  <- apx + (6 * sigma + 12 * tau_right)
  k <- 0L; while (g(lo) > 0 && k < 20L) { lo <- apx - 2 * (apx - lo); k <- k + 1L }
  k <- 0L; while (g(hi) > 0 && k < 20L) { hi <- apx + 2 * (hi - apx); k <- k + 1L }
  L <- tryCatch(stats::uniroot(g, c(lo, apx), tol = 1e-10)$root, error = function(e) NA_real_)
  R <- tryCatch(stats::uniroot(g, c(apx, hi), tol = 1e-10)$root, error = function(e) NA_real_)
  c(apex = apx, fmax = fm, left = L, right = R, width = R - L)
}

#' Fit a sum of BEMG components to one cluster of merged peaks.
#'
#' Variable projection: for given shapes (mu, sigma, tau_l, tau_r) the
#' non-negative component AREAS are solved linearly, so the optimiser only
#' sees 4 nonlinear parameters per component. Parameters are optimised on an
#' unconstrained scale: mu = start + theta * sigma0, sigma and both taus as
#' sigma0 * exp(theta) (clamped to +/- 6, i.e. tau down to 0.25% of sigma --
#' numerically a pure Gaussian side).
#'
#' [MODIFIED] The optimiser is Levenberg-Marquardt on the projected residual,
#' with a finite-difference Jacobian that re-evaluates only the one basis
#' column a parameter touches (k-fold cheaper than a generic optimiser). It is
#' run from `n_restart` starts (tailing, fronting and symmetric shapes; any
#' further starts are jittered) and the best solution is kept; Nelder-Mead +
#' BFGS is the fallback if LM fails. Default 5 starts: the three shapes plus
#' two jittered tailing starts (reproducible under the pipeline's fixed seed).
#'
#' @param starts data.frame(tR, H, sg) initial guesses, one row per component
#'   (tR = apex time, sg = sigma estimate, e.g. FWHM / 2.355).
#' @param tau0 c(left, right) starting time constants in units of sigma.
fit_cluster <- function(
    x,
    y,
    starts,
    tau0      = c(0.05, 0.5),
    maxit     = 100L,
    n_restart = 5L
  ) {
  k   <- nrow(starts)
  tR0 <- starts$tR
  sg0 <- pmax(starts$sg, .Machine$double.eps)
  m   <- length(x)
  cl  <- function(v) pmin(pmax(v, -6), 6)
  par_i <- function(th, i) c(mu = tR0[i] + th[i] * sg0[i],
                             sg = sg0[i] * exp(cl(th[k + i])),
                             tl = sg0[i] * exp(cl(th[2L * k + i])),
                             tr = sg0[i] * exp(cl(th[3L * k + i])))
  col_i <- function(th, i) {
    p <- par_i(th, i)
    .bemg_unit(x, p[["mu"]], p[["sg"]], p[["tl"]], p[["tr"]])
  }
  basis <- function(th) vapply(seq_len(k), function(i) col_i(th, i), numeric(m))

  ## least-squares areas via the k x k normal equations (Cholesky); QR is
  ## the fallback when the Gram matrix is numerically singular
  ls_solve <- function(P) {
    G <- crossprod(P)
    a <- try(as.numeric(backsolve(chol(G), forwardsolve(t(chol(G)), crossprod(P, y)))),
             silent = TRUE)
    if (inherits(a, "try-error") || any(!is.finite(a)))
      a <- try(as.numeric(qr.solve(P, y)), silent = TRUE)
    a
  }
  solveA <- function(Phi) {
    A <- ls_solve(Phi)
    if (inherits(A, "try-error") || any(!is.finite(A))) return(NULL)
    if (any(A < 0)) {
      keep <- A >= 0
      if (!any(keep)) return(rep(0, k))
      A2 <- rep(0, k)
      a  <- ls_solve(Phi[, keep, drop = FALSE])
      if (inherits(a, "try-error") || any(!is.finite(a))) return(NULL)
      A2[keep] <- pmax(a, 0)
      return(A2)
    }
    A
  }
  resid_phi <- function(Phi) {
    if (any(!is.finite(Phi)) || any(colSums(Phi) <= 0)) return(NULL)
    A <- solveA(Phi)
    if (is.null(A)) return(NULL)
    as.numeric(y - Phi %*% A)
  }
  obj <- function(th) {
    r <- resid_phi(basis(th))
    if (is.null(r)) 1e12 else sum(r^2)
  }

  ## Levenberg-Marquardt on the variable-projection residual; `free` selects
  ## the parameters that move (the others stay at their current values)
  lm_vp <- function(th, free = seq_along(th)) {
    Phi <- basis(th); r <- resid_phi(Phi)
    if (is.null(r)) return(NULL)
    f   <- sum(r^2); lam <- 1e-3; p <- length(free); conv <- FALSE
    h   <- 1e-6
    for (it in seq_len(maxit)) {
      J <- matrix(0, m, p)
      for (jj in seq_len(p)) {
        j  <- free[jj]
        i  <- ((j - 1L) %% k) + 1L
        t2 <- th; t2[j] <- t2[j] + h
        P2 <- Phi; P2[, i] <- col_i(t2, i)
        r2 <- resid_phi(P2)
        if (!is.null(r2)) J[, jj] <- (r2 - r) / h
      }
      g   <- crossprod(J, r)
      JtJ <- crossprod(J)
      ## relative floor keeps the damped system well conditioned when a
      ## parameter sits at its box limit (zero Jacobian column)
      dg  <- pmax(diag(JtJ), 1e-8 * max(diag(JtJ), 1e-300))
      improved <- FALSE
      repeat {
        dlt <- try(-solve(JtJ + lam * diag(dg, p), g), silent = TRUE)
        if (inherits(dlt, "try-error") || any(!is.finite(dlt))) {
          lam <- lam * 10; if (lam > 1e12) break; next
        }
        tn <- th; tn[free] <- tn[free] + as.numeric(dlt)
        ## project the log-scale parameters back into the box [-6, 6]
        tn[(k + 1L):(4L * k)] <- cl(tn[(k + 1L):(4L * k)])
        Pn <- basis(tn); rn <- resid_phi(Pn)
        fn <- if (is.null(rn)) Inf else sum(rn^2)
        if (fn < f) {
          rel <- (f - fn) / max(f, .Machine$double.eps)
          th <- tn; Phi <- Pn; r <- rn; f <- fn
          lam <- max(lam / 3, 1e-12); improved <- TRUE
          if (rel < 1e-8 || max(abs(dlt)) < 1e-7) conv <- TRUE
          break
        }
        lam <- lam * 4
        if (lam > 1e12) break
      }
      if (!improved) { conv <- TRUE; break }
      if (conv) break
    }
    list(par = th, value = f, convergence = if (conv) 0L else 1L)
  }

  ## Deterministic starts first: tailing (mu slightly before the apex),
  ## fronting (mirror image) and symmetric; further restarts jitter the
  ## tailing start. A single local optimiser cannot flip the asymmetry of a
  ## peak by itself, so both directions are always tried.
  th_tail <- c(rep(-0.3, k), rep(0, k), rep(log(tau0[1L]), k), rep(log(tau0[2L]), k))
  th_frnt <- c(rep( 0.3, k), rep(0, k), rep(log(tau0[2L]), k), rep(log(tau0[1L]), k))
  th_symm <- c(rep( 0.0, k), rep(0, k), rep(log(0.2), k),     rep(log(0.2), k))
  starts_th <- list(th_tail, th_frnt, th_symm)
  best <- NULL
  for (rr in seq_len(max(1L, n_restart))) {
    th <- if (rr <= length(starts_th)) starts_th[[rr]] else
      th_tail + stats::rnorm(length(th_tail), 0, 0.25)
    ## two stages: first the shapes with every apex held at its seed (stops
    ## components sliding onto each other before their widths are right),
    ## then everything
    f1 <- try(lm_vp(th, free = (k + 1L):(4L * k)), silent = TRUE)
    if (!inherits(f1, "try-error") && !is.null(f1)) th <- f1$par
    f  <- try(lm_vp(th), silent = TRUE)
    if (inherits(f, "try-error") || is.null(f)) {
      f <- stats::optim(th, obj, method = "Nelder-Mead",
                        control = list(maxit = max(800L, 250L * k), reltol = 1e-12))
      f <- stats::optim(f$par, obj, method = "BFGS", control = list(maxit = 800L))
    }
    if (is.null(best) || f$value < best$value) best <- f
  }
  fit <- best

  q    <- lapply(seq_len(k), function(i) par_i(fit$par, i))
  q    <- list(mu = vapply(q, `[[`, 0, "mu"), sg = vapply(q, `[[`, 0, "sg"),
               tl = vapply(q, `[[`, 0, "tl"), tr = vapply(q, `[[`, 0, "tr"))
  Phi  <- basis(fit$par)
  q$A  <- solveA(Phi)
  if (is.null(q$A)) q$A <- rep(0, k)
  dxm  <- mean(diff(x))

  shp <- t(vapply(seq_len(k), function(i)
    bemg_shape(q$mu[i], q$sg[i], q$tl[i], q$tr[i], 0.5), numeric(5)))

  out <- data.frame(
    component   = seq_len(k),
    rt          = shp[, "apex"],
    tR          = shp[, "apex"],
    mu          = q$mu,
    H           = q$A * shp[, "fmax"],
    height      = q$A * shp[, "fmax"],
    fwhm        = shp[, "width"],
    fwhm_left   = shp[, "left"],
    fwhm_right  = shp[, "right"],
    sigma       = q$sg,
    tau_left    = q$tl,
    tau_right   = q$tr,
    area        = q$A,
    area_window = colSums(sweep(Phi, 2L, q$A, "*")) * dxm,
    rss         = fit$value,
    npar        = 5L * k,
    nobs        = length(x),
    converged   = fit$convergence == 0,
    row.names   = NULL
  )
  out <- out[order(out$rt), , drop = FALSE]
  out$component <- seq_len(k)
  rownames(out) <- NULL
  out
}

## --- Main peak finder --------------------------------------------------------
#' Robust peak picking for a baseline-corrected chromatogram
#'
#' HYPERPARAMETERS AND HOW THEY ARE SET (all resolved values are returned in
#' attr(out, "hyperparameters") and copied into the run summary):
#'
#'  noise sd     measured, never tuned: lower quartile of the sd of 201-point
#'               linearly detrended windows (noise_profile()).
#'  fwer         0.05, universal. Probability of >= 1 false peak per run.
#'  prom cutoff  derived, never tuned: (1 - fwer) quantile of the maximum
#'               prominence of `nsim` block-bootstrap noise surrogates.
#'  snr_min      20, universal. Minimum apex height / noise sd. Because sd is
#'               measured per chromatogram, the threshold in intensity units
#'               adapts automatically; the S/N value itself is not optimised
#'               per chromatogram. 20 is deliberately above the ICH LOQ (10)
#'               so that reported retention times and widths are precise.
#'  scales       derived: 12 geometric scales from half the smallest to three
#'               times the largest Gaussian sigma (= FWHM / 2.355, in points)
#'               of the prominence-detected peaks. Only used if do_shoulders.
#'  shoulder_snr derived: largest ridge S/N produced by `nsim_cwt` noise
#'               surrogates run through the identical CWT pipeline. With
#'               nsim_cwt = ceiling(1/fwer) - 1 (19 for fwer = 0.05) the
#'               maximum of the surrogates is an exact Monte-Carlo test at
#'               level fwer.
#'  min_ridge_len 4 of 12 scales: a shoulder must persist over a factor
#'               ~2.3 in scale, which rejects single-scale noise ridges.
#'
#' WHEN ARE SHOULDERS LOOKED FOR? Only when do_shoulders = TRUE. Every CWT
#' ridge is then a candidate; a candidate becomes a "shoulder" only if
#'   (a) it persists over >= min_ridge_len scales,
#'   (b) its coefficient S/N exceeds shoulder_snr,
#'   (c) it is not the ridge of an apex already found by prominence
#'       (distance > max(3, 1.5 x ridge scale) points), and
#'   (d) [NEW] at its own scale the CWT coefficient dips between it and the
#'       curvature lobe of its parent apex by more than shoulder_snr
#'       coefficient-sds, i.e. it is a separate negative-curvature lobe, not
#'       noise riding on the curvature of a single broad apex.
#'
#' RESOLUTION LIMITS. Simulated with this code (two BEMG components of equal
#' width, sigma = 10 points, tau_right = 0.5 sigma, 4 replicates per cell,
#' S/N of the major peak 100 and 1000; Rs from the FWHM). Smallest Rs at which
#' BOTH components were reported in every replicate:
#'
#'   height ratio   prominence   + CWT shoulder   BEMG/BIC decomposition
#'   1 : 1          0.6          0.6              0.4
#'   1 : 0.3        0.8 - 1.0    0.7              0.4 - 0.5
#'   1 : 0.1        1.0 - >1.0   0.8 - 1.0        0.7
#'
#' Nothing was resolved at Rs = 0.3 by any method, and the BIC never proposed
#' a third component. The BEMG column assumes the cluster is fitted, i.e.
#' decomp_min_members = 1 for clusters with a single detected apex; with the
#' default (2) a component hidden under a single detected apex is not looked
#' for. Below Rs ~ 0.35 two components are not identifiable from one channel
#' (select_components(min_fit_rs = 0.35) rejects such fits); only a second
#' data dimension (DAD spectra: peak purity / EFA / MCR-ALS, or MS) can reveal
#' them.
#'
#' @param x,y     retention time and (baseline-corrected) signal.
#' @param fwer    family-wise false-positive rate for the prominence cutoff.
#' @param min_prom  optional hard override of the calibrated cutoff.
#' @param snr_min minimum apex height / noise sd. Applied on top of prominence.
#' @param scales  CWT scales in points. NULL = geometric grid inferred from the
#'                widths of the prominence-detected peaks.
#' @param min_ridge_len  a CWT ridge must persist over this many scales.
#' @param shoulder_snr   a shoulder-only candidate must reach this coefficient
#'                       S/N to be reported. NULL = calibrated on surrogates.
#' @param nsim    bootstrap replicates for the noise calibration.
#' @param nsim_cwt surrogates for the shoulder threshold. NULL =
#'                ceiling(1 / fwer) - 1.
#' @param denoise [NEW] denoiser applied after baseline correction, one of
#'                DENOISE_METHODS. Default "none".
#' @param denoise_args [NEW] list passed to denoise_chromatogram(): param,
#'                tune ("rv" | "acf" | "none"), grid, sg_order.
#' @param noise_win,noise_q [NEW, exposed] noise_profile() window (points)
#'                and quantile.
#' @param shoulder_valley [NEW] apply test (d) above.
#'
#' @return data.frame, one row per peak:
#'   rt, idx, height, prominence, snr, fwhm, left/right (integration bounds),
#'   area, type ("peak" | "shoulder"), cluster, resolution (Rs to the previous
#'   peak in the same cluster). attr "denoised" holds the detector signal when
#'   a denoiser was used.
pick_peaks <- function(
    x,
    y,
    fwer            = 0.05,
    min_prom        = NULL,
    snr_min         = 20,
    scales          = NULL,
    min_ridge_len   = 4L,
    shoulder_snr    = NULL,
    nsim            = 25L,
    nsim_cwt        = NULL,
    do_shoulders    = TRUE,
    denoise         = "none",
    denoise_args    = list(),
    noise_win       = 201L,
    noise_q         = 0.25,
    shoulder_valley = TRUE
  ) {

  n     <- length(y)
  dx    <- mean(diff(x))
  y_raw <- y

  ## noise of the RAW (baseline-corrected) trace: defines S/N, the snr_min
  ## floor and the surrogate pool, exactly as a chromatographer measures it
  np    <- noise_profile(y_raw, win = noise_win, q = noise_q)
  sd_n  <- np$sd

  ## ---- (0) optional denoising -----------------------------------------------
  denoise <- match.arg(denoise, DENOISE_METHODS)
  dn      <- NULL
  den_fun <- NULL
  if (denoise != "none") {
    dn <- do.call(denoise_chromatogram,
                  c(list(y = y_raw, method = denoise), denoise_args))
    y       <- dn$y
    den_fun <- .dn_make_applier(dn$method, dn$param, dn$sg_order)
  }

  ## ---- (1)+(2) prominence, with a bootstrapped cutoff ---------------------
  cal  <- calibrate_prominence(y_raw, nsim = nsim, fwer = fwer, np = np,
                               denoiser = den_fun)
  if (!is.null(min_prom)) cal$cutoff <- min_prom
  pr   <- prominence(y)
  keep <- pr$prominence >= cal$cutoff & y[pr$idx] >= snr_min * sd_n
  pk   <- pr[keep, , drop = FALSE]

  ## ---- widths, needed to choose CWT scales --------------------------------
  fwhm_crossings <- function(i, lo, hi) {
    h <- y[i] / 2
    l <- i; while (l > lo && y[l] > h) l <- l - 1L
    r <- i; while (r < hi && y[r] > h) r <- r + 1L

    list(
      fwhm       = (r - l) * dx,
      fwhm_left  = l,
      fwhm_right = r
    )
  }
  if (nrow(pk)) {
    fwhm_result   <- mapply(fwhm_crossings, pk$idx, pk$left, pk$right)
    pk$fwhm       <- as.numeric(fwhm_result["fwhm", ])
    pk$fwhm_left  <- as.numeric(fwhm_result["fwhm_left", ])
    pk$fwhm_right <- as.numeric(fwhm_result["fwhm_right", ])
  } else {
    pk$fwhm_left <- pk$fwhm_right <- pk$fwhm <- numeric(0)
  }

  ## sigma (points) implied by each FWHM, for a Gaussian
  sig_pts <- if (nrow(pk)) pmax(1, pk$fwhm / dx / 2.355) else 3
  scales_auto <- is.null(scales)
  if (scales_auto) {
    lo <- max(1, min(sig_pts) / 2); hi <- max(sig_pts) * 3
    scales <- unique(round(exp(seq(log(lo), log(max(hi, lo * 4)), length.out = 12)), 2))
  }

  ## ---- (3) CWT ridges ------------------------------------------------------
  sh <- NULL
  if (is.null(nsim_cwt)) nsim_cwt <- max(5L, as.integer(ceiling(1 / fwer)) - 1L)
  shoulder_snr_used <- NA_real_
  if (do_shoulders) {
    cw <- cwt_ridges(y, scales)
    if (!is.null(cw$ridges)) {
      rid <- cw$ridges[cw$ridges$len >= min_ridge_len, , drop = FALSE]
      ## Noise level of each coefficient, computed exactly from the empirical
      ## noise autocovariance (after denoising) rather than assuming white noise.
      csd_of <- function(s) .coef_sd(.ricker(s), cal$gamma)
      rid$csd <- vapply(rid$scale, csd_of, numeric(1))
      rid$snr <- rid$coef / rid$csd
      ## Threshold calibrated the same way as the prominence cutoff: run the
      ## noise surrogates through the identical (denoise +) CWT + ridge
      ## pipeline and take the largest ridge S/N noise produces. Any bias in
      ## .coef_sd cancels, and it accounts for the multiple testing across
      ## ~n/(2s) positions and every scale.
      if (is.null(shoulder_snr)) {
        b  <- cal$block
        mx <- vapply(seq_len(nsim_cwt), function(i) {
          sg <- rep_len(.surrogate(np$pool, b), n)
          if (!is.null(den_fun)) sg <- den_fun(sg)
          cs <- cwt_ridges(sg, scales)
          if (is.null(cs$ridges)) return(0)
          rr <- cs$ridges[cs$ridges$len >= min_ridge_len, , drop = FALSE]
          if (!nrow(rr)) return(0)
          max(rr$coef / vapply(rr$scale, csd_of, numeric(1)))
        }, numeric(1))
        shoulder_snr <- max(mx)
      }
      shoulder_snr_used <- shoulder_snr
      rid <- rid[which(rid$snr >= shoulder_snr), , drop = FALSE]
      ## drop ridges that merely re-detect an already-found apex
      if (nrow(rid) && nrow(pk)) {
        tol  <- pmax(3, 1.5 * rid$scale)
        dmin <- vapply(seq_len(nrow(rid)),
                       function(j) min(abs(rid$idx_fine[j] - pk$idx)), numeric(1))
        rid  <- rid[dmin > tol, , drop = FALSE]
      }
      ## [NEW] curvature-valley test. At the ridge's own scale the CWT
      ## coefficient is a smoothed -y''. A genuine shoulder is a SEPARATE
      ## negative-curvature lobe: between it and the lobe of its parent apex
      ## the coefficient must dip by more than shoulder_snr coefficient-sds.
      ## The parent lobe is the largest coefficient within the parent's
      ## half-width (excluding the half of the gap nearest the candidate).
      ## Noise riding on the curvature hump of one broad apex -- which can
      ## have a large ridge S/N because the hump itself is large -- never
      ## produces such a dip and is rejected.
      if (nrow(rid) && nrow(pk) && isTRUE(shoulder_valley)) {
        js  <- match(rid$scale, cw$scales)
        okv <- logical(nrow(rid))
        for (r in seq_len(nrow(rid))) {
          s_i <- rid$idx[r]
          a   <- which.min(abs(pk$idx - s_i))
          a_i <- pk$idx[a]
          hwa <- max(2L, as.integer(ceiling(pk$fwhm[a] / dx / 2)))
          win <- max(1L, a_i - hwa):min(n, a_i + hwa)
          win <- win[abs(win - s_i) > abs(a_i - s_i) / 2]
          if (!length(win)) { okv[r] <- TRUE; next }
          wj  <- cw$W[, js[r]]
          p_i <- win[which.max(wj[win])]
          seg <- wj[min(s_i, p_i):max(s_i, p_i)]
          dip <- min(wj[s_i], wj[p_i]) - min(seg)
          okv[r] <- is.finite(dip) && dip >= shoulder_snr * rid$csd[r]
        }
        rid <- rid[okv, , drop = FALSE]
      }
      if (nrow(rid)) sh <- rid
    }
  }

  ## ---- assemble ------------------------------------------------------------
  out <- data.frame(
    idx        = integer(0),
    height     = numeric(0),
    prominence = numeric(0),
    fwhm       = numeric(0),
    fwhm_left  = numeric(0),
    fwhm_right = numeric(0),
    type       = character(0)
  )
  if (nrow(pk))
    out <- rbind(out,
                 data.frame(
                   idx        = pk$idx,
                   height     = y[pk$idx],
                   prominence = pk$prominence,
                   fwhm       = pk$fwhm,
                   fwhm_left  = x[pk$fwhm_left],
                   fwhm_right = x[pk$fwhm_right],
                   type       = "peak")
                 )
  if (!is.null(sh) && nrow(sh))
    out <- rbind(out,
                 data.frame(
                   idx        = sh$idx_fine,
                   height     = y[sh$idx_fine],
                   prominence = NA_real_,
                   fwhm       = sh$scale * 2.355 * dx,
                   fwhm_left  = NA_real_,
                   fwhm_right = NA_real_,
                   type       = "shoulder")
                 )

  hyper <- list(
    noise_win        = np$win,
    noise_q          = np$q,
    noise_sd         = sd_n,
    noise_n_quiet    = np$n_quiet,
    noise_acf_len    = cal$acf_len,
    bootstrap_block  = cal$block,
    fwer             = fwer,
    nsim             = nsim,
    prom_cutoff      = cal$cutoff,
    prom_cutoff_auto = is.null(min_prom),
    snr_min          = snr_min,
    height_floor     = snr_min * sd_n,
    do_shoulders     = do_shoulders,
    scales           = scales,
    scales_auto      = scales_auto,
    min_ridge_len    = min_ridge_len,
    nsim_cwt         = if (do_shoulders) nsim_cwt else NA_integer_,
    shoulder_snr     = shoulder_snr_used,
    shoulder_valley  = shoulder_valley,
    denoise          = if (is.null(dn)) "none" else dn$method,
    denoise_param    = if (is.null(dn)) NA_real_ else dn$param,
    denoise_tune     = if (is.null(dn)) NA_character_ else dn$tune,
    denoise_target   = if (is.null(dn)) NA_real_ else dn$target,
    denoise_kernel_fwhm_pts = if (is.null(dn)) 0 else dn$kernel_fwhm_pts,
    denoise_kernel_fwhm_min = if (is.null(dn)) 0 else dn$kernel_fwhm_pts * dx
  )

  set_attrs <- function(o) {
    attr(o, "noise_sd")        <- sd_n
    attr(o, "prom_cutoff")     <- cal$cutoff
    attr(o, "scales")          <- scales
    attr(o, "null_max")        <- cal$null_max
    attr(o, "shoulder_snr")    <- if (do_shoulders) shoulder_snr_used else NA
    attr(o, "denoised")        <- if (is.null(dn)) NULL else y
    attr(o, "hyperparameters") <- hyper
    o
  }

  if (!nrow(out)) return(set_attrs(out))
  out <- out[order(out$idx), ]
  ## edge guard: prominence is ill-defined against a truncated flank
  edge <- max(3L, ceiling(min(scales)))
  out  <- out[out$idx > edge & out$idx <= n - edge, , drop = FALSE]
  if (!nrow(out)) return(set_attrs(out))

  ## integration bounds: valley between neighbours, or the point where the
  ## signal falls into the noise, whichever comes first
  m   <- nrow(out)
  bnd <- matrix(NA_integer_, m, 2)
  for (j in seq_len(m)) {
    i  <- out$idx[j]
    lo <- if (j > 1L) out$idx[j - 1L] else 1L
    hi <- if (j < m)  out$idx[j + 1L] else n
    vl <- if (j > 1L) lo + which.min(y[lo:i]) - 1L else 1L
    vr <- if (j < m)  i  + which.min(y[i:hi]) - 1L else n
    ## walk out from the apex; stop at the valley or when the trace drops
    ## into the noise, whichever comes first
    l <- i; while (l > vl && y[l] > snr_min * sd_n) l <- l - 1L
    r <- i; while (r < vr && y[r] > snr_min * sd_n) r <- r + 1L
    bnd[j, ] <- c(l, r)
  }
  out$left <- bnd[, 1]; out$right <- bnd[, 2]

  ## One more security check of FWHM if bands moved
  if (nrow(out)) {
    fwhm_result   <- mapply(fwhm_crossings, out$idx, out$left, out$right)
    out$fwhm       <- as.numeric(fwhm_result["fwhm", ])
    out$fwhm_left  <- x[as.numeric(fwhm_result["fwhm_left", ])]
    out$fwhm_right <- x[as.numeric(fwhm_result["fwhm_right", ])]
  }

  ## [MODIFIED] Saturated (detector-clipped) apexes. A clipped peak is flat
  ## in the RAW data to well within the noise: >= 4 consecutive points within
  ## 0.1 noise-sd of the apex value. The old rule (points >= 99.5% of the apex)
  ## also fired on the genuinely flat top of two merged, equal peaks, which
  ## then blocked their decomposition. Shoulders found inside a clipped
  ## plateau are curvature artefacts of the plateau edges and are dropped.
  plateau_tol <- max(0.1 * sd_n, 1e-9 * max(abs(y_raw)))
  out$saturated <- vapply(seq_len(m), function(j) {
    lo_j <- out$left[j]; hi_j <- out$right[j]
    seg  <- y_raw[lo_j:hi_j]
    top  <- max(seg)
    r    <- rle(abs(seg - top) <= plateau_tol)
    any(r$values & r$lengths >= 4L)
  }, logical(1))
  if (any(out$saturated) && any(out$type == "shoulder")) {
    sat <- which(out$saturated & out$type == "peak")
    bad <- vapply(seq_len(m), function(j) out$type[j] == "shoulder" &&
                    any(out$idx[j] >= out$left[sat] & out$idx[j] <= out$right[sat]),
                  logical(1))
    if (any(bad)) {
      out <- out[!bad, , drop = FALSE]; m <- nrow(out)
      bnd <- cbind(out$left, out$right)
    }
  }
  out$rt   <- x[out$idx]
  out$snr  <- out$height / sd_n
  out$area <- vapply(seq_len(m), function(j)
    sum(y[out$left[j]:out$right[j]]) * dx, numeric(1))

  ## cluster: peaks whose integration ranges touch are one merged group
  cl <- integer(m); cl[1] <- 1L
  if (m > 1L) for (j in 2:m)
    cl[j] <- if (out$left[j] <= out$right[j - 1L]) cl[j - 1L] else cl[j - 1L] + 1L
  out$cluster <- cl

  ## chromatographic resolution to the previous peak in the same cluster
  out$resolution <- NA_real_
  if (m > 1L) for (j in 2:m) if (cl[j] == cl[j - 1L]) {
    w <- out$fwhm[j] + out$fwhm[j - 1L]
    out$resolution[j] <- if (w > 0) 1.18 * (out$rt[j] - out$rt[j - 1L]) / w else NA
  }

  rownames(out) <- NULL
  out <- out[, c("rt",   "idx",       "height",     "prominence", "snr",
                 "fwhm", "fwhm_left", "fwhm_right", "left",       "right",
                 "area", "type",      "saturated",  "cluster",    "resolution")]
  set_attrs(out)
}

#' How many components does this cluster actually need?
#'
#' Fits k = 1..k_max BEMG profiles and selects k by BIC. Extra components are
#' seeded at the largest positive residual of the previous fit, so a merged
#' peak that is invisible in both the trace AND its second derivative can still
#' be recovered -- which matters, because a shoulder at Rs ~ 0.45 with a 3:1
#' height ratio produces no separate maximum in y and none in y'' either.
#'
#' The BIC uses an EFFECTIVE sample size m / L, where L is the noise
#' autocorrelation length. Chromatographic noise is not white; using the raw
#' point count would make BIC think it has far more evidence than it does and
#' it would happily fit components to wiggles. Each BEMG component costs 5
#' parameters (area, mu, sigma, tau_left, tau_right).
#'
#' @param model_tol [NEW] relative model-error floor (fraction of the cluster
#'   maximum) added to the residual variance in the BIC.
#' @param min_fit_rs [NEW] models in which two fitted components are closer
#'   than this resolution are rejected as unidentifiable.
select_components <- function(
    x,
    y,
    seeds,
    k_max         = 3L,
    acl           = 1,
    bic_margin    = 10,
    min_area_frac = 0.01,
    model_tol     = 0.01,
    min_fit_rs    = 0.35
  ) {
  m     <- length(y)
  m_eff <- max(5, m / max(1, acl))
  ## [NEW] model-error floor: no analytical peak shape describes a real peak
  ## to better than ~1% of its height. Without this floor, at S/N >~ 500 the
  ## BIC treats every percent of shape mismatch as evidence for another
  ## component and splits single peaks. The floor enters as an extra variance
  ## (model_tol * max(y))^2 added to rss / m.
  v_mod <- (model_tol * max(y))^2
  best  <- NULL; best_bic <- Inf; hist <- NULL
  st    <- seeds[order(-seeds$H), , drop = FALSE]
  prev  <- NULL
  for (k in seq_len(k_max)) {
    if (k <= nrow(st)) {
      f <- try(fit_cluster(x, y, st[seq_len(k), , drop = FALSE]), silent = TRUE)
    } else {
      ## Two ways to propose one more component, both tried, lower RSS kept:
      ##  (a) add a component at the largest positive residual of the last
      ##      fit (a hidden peak beside a larger one);
      ##  (b) split the widest component of the best admissible fit into two
      ##      at +/- FWHM/4 (two near-equal peaks forming one flat top, which
      ##      leaves no single residual maximum).
      cand <- list()
      base <- if (!is.null(prev)) prev else best
      if (!is.null(base)) {
        pred <- numeric(m)
        for (i in seq_len(nrow(base)))
          pred <- pred + bemg(x, base$area[i], base$mu[i], base$sigma[i],
                              base$tau_left[i], base$tau_right[i])
        res <- y - pred
        j <- which.max(res)
        if (res[j] > 0) cand[[length(cand) + 1L]] <- rbind(
          data.frame(tR = base$rt, H = base$height, sg = base$sigma),
          data.frame(tR = x[j],    H = res[j],      sg = stats::median(base$sigma)))
      }
      bb <- if (!is.null(best)) best else base
      if (!is.null(bb) && nrow(bb) == k - 1L) {
        w  <- which.max(bb$fwhm)
        sp <- data.frame(tR = bb$rt[w] + c(-1, 1) * bb$fwhm[w] / 4,
                         H  = bb$height[w], sg = bb$sigma[w] * 0.6)
        cand[[length(cand) + 1L]] <- rbind(
          data.frame(tR = bb$rt[-w], H = bb$height[-w], sg = bb$sigma[-w]), sp)
      }
      if (!length(cand)) break
      fits <- lapply(cand, function(s0) try(fit_cluster(x, y, s0), silent = TRUE))
      fits <- fits[!vapply(fits, inherits, logical(1), "try-error")]
      if (!length(fits)) next
      f <- fits[[which.min(vapply(fits, function(z) z$rss[1], numeric(1)))]]
    }
    if (inherits(f, "try-error")) next
    prev <- f
    kk   <- nrow(f)
    rs_min <- if (kk > 1L) {
      min(1.18 * diff(f$rt) / (f$fwhm[-1L] + f$fwhm[-kk]), na.rm = TRUE)
    } else NA_real_
    bic  <- m_eff * log(f$rss[1] / m + v_mod) + 5 * kk * log(m_eff)
    ## a component carrying a negligible share of the cluster, or two
    ## components closer than min_fit_rs (not identifiable from one channel,
    ## the signature of the optimiser bending two shapes into one peak), make
    ## the whole model inadmissible rather than reportable
    a_frac <- if (sum(f$area) > 0) min(f$area) / sum(f$area) else NA_real_
    admissible <- isTRUE(all(is.finite(f$rt))) && isTRUE(sum(f$area) > 0) &&
      !(kk > 1L && (!is.finite(a_frac) || a_frac < min_area_frac ||
                    !is.finite(rs_min) || rs_min < min_fit_rs))
    hist <- rbind(hist, data.frame(k = kk, rss = f$rss[1], bic = bic,
                                   min_rs = rs_min, admissible = admissible))
    if (!admissible) next
    ## accept a richer model only on strong evidence (Kass & Raftery: a BIC
    ## drop of >10 is "very strong"). Without this the fit happily adds a
    ## component to every bit of tailing it cannot otherwise explain.
    if (bic < best_bic - (if (is.null(best)) 0 else bic_margin)) {
      best_bic <- bic; best <- f
    }
  }
  list(fit = best, bic = hist)
}

#' Refit every multi-component cluster with a sum of BEMG profiles.
#' Returns one row per fitted component with apex, height, FWHM and area.
#'
#' @param k_max  minimum upper bound on the number of components tried.
#' @param k_extra [NEW] the search always goes up to (detected members +
#'   k_extra) components, so up to k_extra hidden components per cluster can
#'   be proposed by BIC (was fixed at +1).
decompose_clusters <- function(
    x,
    y,
    pk,
    min_members = 1L,
    pad         = 1.5,
    k_max       = 3L,
    k_extra     = 2L,
    min_snr     = 20,
    acl         = NULL,
    model_tol   = 0.01,
    min_fit_rs  = 0.35
  ) {
  if (is.null(acl)) {
    np <- noise_profile(y); acl <- .acf_len(np$pool)
  }
  dxm <- mean(diff(x))
  res <- list(); bics <- list()
  for (cid in unique(pk$cluster)) {
    sub <- pk[pk$cluster == cid, , drop = FALSE]
    if (nrow(sub) < min_members) next
    if (max(sub$snr) < min_snr)  next             # not worth fitting
    if (any(sub$saturated))      next             # a clipped peak has no shape
    lo  <- max(1L, min(sub$left)  - round(pad * mean(sub$fwhm) / dxm))
    hi  <- min(length(x), max(sub$right) + round(pad * mean(sub$fwhm) / dxm))
    st  <- data.frame(tR = sub$rt, H = sub$height, sg = pmax(sub$fwhm, 1e-6) / 2.355)
    sel <- try(select_components(x[lo:hi], y[lo:hi], st,
                                 k_max = max(k_max, nrow(st) + k_extra), acl = acl,
                                 model_tol = model_tol, min_fit_rs = min_fit_rs),
               silent = TRUE)
    if (inherits(sel, "try-error") || is.null(sel$fit)) next
    f <- sel$fit; f$cluster <- cid
    f$k_detected <- nrow(st); f$k_fitted <- nrow(f)
    res[[length(res) + 1L]] <- f
    sel$bic$cluster <- cid; bics[[length(bics) + 1L]] <- sel$bic
  }
  if (!length(res)) return(NULL)
  out <- do.call(rbind, res)
  attr(out, "bic") <- do.call(rbind, bics)
  out
}

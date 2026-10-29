## ---------------------------------------------------------------------------
##  FlatFit -- Flatness-weighted baseline fit.
##
##  Port of, and extension to, the FlatFit algorithm from MOCCA2:
##
##    J. Obořil, C. P. Haas, M. Lübbesmeyer, R. Nicholls, T. Gressling,
##    K. F. Jensen, G. Volpin, J. Hillenbrand,
##    "Automated processing of chromatograms: a comprehensive python package
##     with a GUI for intelligent peak identification and deconvolution in
##     chemical reaction analysis", Digital Discovery 3 (2024) 2041-2051.
##    doi:10.1039/D4DD00214H
##
##  While AsLS and arPLS decide points belonging to the baseline from the residual
##  y - b, so the weights depend on the current baseline and the algorithm has
##  to be iterated. FlatFit instead decides from the FLATNESS of the data --
##  small first and second derivative -- which does not depend on b at all. The
##  weights can therefore be computed once and the penalised least squares
##  problem
##
##      minimise  (y - b)' W (y - b)  +  lambda * b' D2' D2 b
##
##  has a single closed-form solution without the need for iteration:
##                      (W + lambda D2'D2) b = W y
##
##  [MODIFIED v0.6.0] This file now also holds the iterative members of the
##  same Whittaker family -- PSALSA, arPLS and AsLS (whittaker_baseline(),
##  bottom of the file) -- ported from the pybaselines definitions, with an
##  automatic, sampling-rate-invariant choice of lambda from the widest peak.
##  FlatFit's lambda: lambda = smoothness * n^4 (MOCCA2 scaling, smoothness =
##  1 by default, a preset, not fitted); its Savitzky-Golay half-window is
##  estimated from the data (.ff_auto_window(), see below).
## ---------------------------------------------------------------------------

# ---------- Small helpers  ----------------------------------------------------

.ff_pad <- function(y, hw) {
  n <- length(y)
  if (hw < 1) return(y)
  hw <- min(hw, n - 1L)
  c(rev(y[2:(hw + 1L)]), y, rev(y[(n - hw):(n - 1L)]))
}

.ff_conv <- function(y, k) {
  hw <- (length(k) - 1L) %/% 2L
  if (hw < 1L) return(y * sum(k))
  n <- length(y); yp <- .ff_pad(y, hw)
  out <- numeric(n)
  for (j in seq_along(k)) out <- out + k[j] * yp[j:(j + n - 1L)]
  out
}

## Savitzky-Golay coefficients, polynomial order `ord`, derivative `d`.
## Note the scaled basis u = b / hw. Built on raw powers of b, the Vandermonde
## becomes singular for long windows -- MOCCA2's default p = 0.03 on a
## 24 000-point trace asks for a 721-point window, where b^3 reaches 5e7 and
## crossprod(X) has reciprocal condition number 3e-18. Working in u and
## rescaling the derivative by hw^-d is exact and unconditionally stable.
.ff_sg_coef <- function(hw, ord = 3L, d = 0L) {
  u <- (-hw:hw) / hw
  X <- outer(u, 0:ord, "^")
  C <- solve(crossprod(X), t(X))
  factorial(d) * C[d + 1L, ] / hw^d
}

## Savitzky-Golay filter with scipy's mode = "interp" edge handling: instead of
## padding, fit ONE polynomial of the given order to the first (and last) full
## window and evaluate it at the edge positions.
.ff_sg <- function(y, hw, ord = 3L, d = 0L, edge = c("interp", "reflect")) {
  edge <- match.arg(edge)
  if (hw < ord) hw <- ord + 1L
  n <- length(y)
  win <- 2L * hw + 1L
  if (edge == "reflect" || n < win) return(.ff_conv(y, .ff_sg_coef(hw, ord, d)))

  out <- .ff_conv(y, .ff_sg_coef(hw, ord, d))

  uu <- ((0:(win - 1L)) - hw) / hw          # scaled, as above
  X  <- outer(uu, 0:ord, "^")
  XtXi_Xt <- solve(crossprod(X), t(X))
  ## derivative of sum_k beta_k u^k w.r.t. the ORIGINAL variable
  dcoef <- function(beta, u0) {
    if (d == 0L) return(as.numeric(outer(u0, 0:ord, "^") %*% beta))
    k <- d:ord
    if (!length(k)) return(rep(0, length(u0)))
    fac <- vapply(k, function(kk) prod((kk - d + 1L):kk), numeric(1))
    as.numeric(outer(u0, k - d, "^") %*% (fac * beta[k + 1L])) / hw^d
  }
  bL <- as.numeric(XtXi_Xt %*% y[seq_len(win)])
  out[seq_len(hw)] <- dcoef(bL, uu[seq_len(hw)])
  bR <- as.numeric(XtXi_Xt %*% y[(n - win + 1L):n])
  out[(n - hw + 1L):n] <- dcoef(bR, uu[(hw + 2L):win])
  out
}

## numpy.gradient: central differences inside, one-sided at the ends
.ff_gradient <- function(v) {
  n <- length(v)
  if (n < 2L) return(rep(0, n))
  g <- numeric(n)
  g[1] <- v[2] - v[1]
  g[n] <- v[n] - v[n - 1L]
  if (n > 2L) g[2:(n - 1L)] <- (v[3:n] - v[1:(n - 2L)]) / 2
  g
}

## Symmetric pentadiagonal Cholesky solve, O(n). a0 main diagonal (n),
## a1 first sub-diagonal (n-1), a2 second (n-2).
.ff_penta_solve <- function(a0, a1, a2, b) {
  n <- length(a0)
  d <- numeric(n); e <- numeric(n); f <- numeric(n)
  for (i in seq_len(n)) {
    em <- if (i >= 2L) e[i - 1L] else 0
    fm <- if (i >= 3L) f[i - 2L] else 0
    v <- a0[i] - em * em - fm * fm
    if (!is.finite(v) || v <= 0) v <- .Machine$double.eps
    d[i] <- sqrt(v)
    fm1 <- if (i >= 2L) f[i - 1L] else 0
    if (i <= n - 1L) e[i] <- (a1[i] - fm1 * em) / d[i]
    if (i <= n - 2L) f[i] <- a2[i] / d[i]
  }
  z <- numeric(n)
  for (i in seq_len(n)) {
    s <- b[i]
    if (i >= 2L) s <- s - e[i - 1L] * z[i - 1L]
    if (i >= 3L) s <- s - f[i - 2L] * z[i - 2L]
    z[i] <- s / d[i]
  }
  x <- numeric(n)
  for (i in rev(seq_len(n))) {
    s <- z[i]
    if (i <= n - 1L) s <- s - e[i] * x[i + 1L]
    if (i <= n - 2L) s <- s - f[i] * x[i + 2L]
    x[i] <- s / d[i]
  }
  x
}

## Whittaker solve: (W + lambda * D2'D2) b = W y, W diagonal.
## t(D2) %*% D2 for the second-difference operator has bands
##   diag  1, 5, 6, 6, ..., 6, 5, 1
##   off1 -2, -4, ..., -4, -2
##   off2  1, ..., 1
## Conditioning of (W + lambda D2'D2) is governed by lambda / min(w). MOCCA2
## scales lambda = smoothness * n^4.
.ff_check_conditioning <- function(w, lambda, what = "baseline fit") {
  r <- lambda / max(stats::median(w), .Machine$double.eps)
  if (r > 1e15) {
    warning(sprintf(
      "%s: lambda / median(w) = %.1e exceeds what double precision can resolve.\n  Reduce `smoothness` (try %.1g) or use fewer points.",
      what, r, 1e15 / r), call. = FALSE)
  }
  invisible(r)
}

.ff_whittaker <- function(y, w, lambda) {
  n <- length(y)
  if (n < 5L) return(rep(stats::median(y), n))
  g0 <- c(1, 5, rep(6, n - 4L), 5, 1)
  g1 <- c(-2, rep(-4, n - 3L), -2)
  g2 <- rep(1, n - 2L)
  .ff_penta_solve(w + lambda * g0, lambda * g1, lambda * g2, w * y)
}

# ---------- The weights -------------------------------------------------------

#' Flatness weights.
#'
#' The MOCCA2 form is
#'
#'   slope     = SG(y, window, order 3, deriv 1)
#'   curvature = gradient(slope)
#'   s = slope^2     / sum(slope^2)         # normalised, sums to 1
#'   c = curvature^2 / sum(curvature^2)
#'   w = 1 / (s + c + eps)
#'
#' Both terms are normalised to unit total, so the two derivative orders
#' contribute comparably regardless of their physical units, and `eps` caps the
#' weight of a perfectly flat point at 1/eps rather than infinity. Nothing here
#' depends on the baseline, which is what makes the fit closed-form.
#'
#' @param eps floor preventing division by zero (MOCCA2 uses 1e-10; the paper's
#'   Table 1 quotes 1e-7 -- the code is authoritative, and the value only sets
#'   the ceiling on how strongly a perfectly flat point can be anchored).
flatfit_weights <- function(y, sg_hw, sg_order = 3L, eps = 1e-10) {
  slope <- .ff_sg(y, sg_hw, ord = sg_order, d = 1L)
  curv  <- .ff_gradient(slope)
  s <- slope^2; ss <- sum(s);  if (ss > 0) s <- s / ss
  c2 <- curv^2; sc <- sum(c2); if (sc > 0) c2 <- c2 / sc
  1 / (s + c2 + eps)
}

# ---------- Main --------------------------------------------------------------

#' Flatness-weighted baseline estimation
#'
#' @param y numeric vector, the measured signal.
#' @param x optional abscissa. Used only for the returned fit and for the
#'   uniform-spacing check; the penalty is in sample units, as in MOCCA2.
#' @param smoothness penalty size. lambda = smoothness * n^4, the MOCCA2
#'   scaling, which makes the result approximately invariant to sampling rate
#'   at fixed run length.
#' @param p Savitzky-Golay window as a fraction of the trace length, MOCCA2
#'   style (default 0.03). Ignored when `sg_hw` is given. See the note below.
#' @param sg_hw Savitzky-Golay HALF-window in points. Overrides `p`. Setting
#'   this to roughly one peak FWHM is far more robust than `p` on long,
#'   densely-sampled traces.
#' @param sg_from_peak_width if TRUE (default) and `sg_hw` is NULL, estimate the
#'   window from the data's own feature scale instead of from `p * n`.
#' @param direction "up" for peaks above the baseline, "down" for absorption
#'   dips, "both" to leave the fit sign-blind (native FlatFit).
#' @param refine number of asymmetric reweighting passes after the closed-form
#'   solve, to remove the documented tendency to sit too high in crowded
#'   regions. 0 reproduces MOCCA2 exactly.
#' @param variant "mocca" for a line-for-line port (forces
#'   sg_from_peak_width = FALSE, refine = 0, direction = "both"), or "extended".
#'
#' @return list with `baseline`, `corrected`, `weights`, `lambda`, `sg_hw`,
#'   `iterations`.
flatfit <- function(
    y,
    x                  = seq_along(y),
    smoothness         = 1,
    p                  = 0.03,
    sg_hw              = NULL,
    sg_from_peak_width = TRUE,
    sg_order           = 3L,
    eps                = 1e-10,
    direction          = c("up", "down", "both"),
    refine             = 3L,
    refine_tol         = 1e-3,
    clamp              = FALSE,
    variant            = c("extended", "mocca")
  ) {

  direction <- match.arg(direction)
  variant   <- match.arg(variant)
  stopifnot(length(x) == length(y), all(is.finite(y)))
  n <- length(y)
  if (n < 5L) stop("flatfit() needs at least 5 points.")

  if (variant == "mocca") {
    sg_from_peak_width <- FALSE
    refine <- 0L
    direction <- "both"
  }

  ## ---- Savitzky-Golay window ----------------------------------------------
  if (is.null(sg_hw)) {
    if (sg_from_peak_width) {
      sg_hw <- .ff_auto_window(y)
    } else {
      ## MOCCA2: filter_window = max(4, int(L * p)); make it odd and take the
      ## half-width. scipy needs an odd window; the reference can hand it an
      ## even one, which errors on current scipy -- rounded up here.
      fw <- max(4L, as.integer(n * p))
      if (fw %% 2L == 0L) fw <- fw + 1L
      sg_hw <- (fw - 1L) %/% 2L
    }
  }
  sg_hw <- max(as.integer(sg_hw), sg_order + 1L)

  ## ---- weights and the closed-form solve ----------------------------------
  w <- flatfit_weights(y, sg_hw = sg_hw, sg_order = sg_order, eps = eps)
  lambda <- smoothness * n^4

  .ff_check_conditioning(w, lambda, "flatfit()")
  z <- .ff_whittaker(y, w, lambda)

  ## ---- optional asymmetric refinement -------------------------------------
  ##
  ## The paper's own stated weakness: "FlatFit overestimates the baseline in
  ## areas with many peaks". The mechanism is visible in the weights -- a
  ## valley between two overlapping peaks has slope ~ 0 and only moderate
  ## curvature, so it scores as "flat" and gets anchored as if it were
  ## baseline. Flatness alone cannot tell a valley floor from a real baseline.
  ##
  ## The residual can. After the closed-form solve, multiply the flatness
  ## weights by an arPLS-style sigmoid that discounts points sitting on the
  ## peak side of the current estimate. This keeps FlatFit's virtues -- the
  ## first pass is still closed-form and parameter-free, so the refinement
  ## starts from a good baseline and converges in 2-3 passes instead of the
  ## dozens arPLS needs from a cold start.
  it <- 0L
  if (refine > 0L && direction != "both") {
    sgn <- if (direction == "up") 1 else -1
    rng <- max(abs(z)); if (!is.finite(rng) || rng <= 0) rng <- 1
    for (k in seq_len(refine)) {
      d <- sgn * (y - z)
      neg <- d[d < 0]
      if (length(neg) < 3L) break
      m <- mean(neg); s <- stats::sd(neg)
      if (!is.finite(s) || s <= 0) break
      wa <- 1 / (1 + exp(2 * (d - (2 * s - m)) / s))
      z_new <- .ff_whittaker(y, w * wa, lambda)
      chg <- max(abs(z_new - z)) / rng
      z <- z_new
      it <- k
      if (chg < refine_tol) break
    }
  }

  if (clamp) z <- if (direction == "down") pmax(z, y) else pmin(z, y)

  list(
    baseline   = z,
    corrected  = y - z,
    weights    = w,
    lambda     = lambda,
    smoothness = smoothness,
    sg_hw      = sg_hw,
    iterations = it,
    variant    = variant
  )
}

# ---------- Automatic Savitzky-Golay window -----------------------------------

## MOCCA2 sets the derivative window to p * n, i.e. to a fraction of the TRACE
## LENGTH. That is safe for the short, modestly-sampled runs it was built for
## (n ~ 1-3k), where 0.03n lands near one peak width by coincidence. On a
## 20-minute run at 20 Hz, 0.03n is 720 points -- twenty times a peak width --
## and a Savitzky-Golay derivative over that window flattens the peaks
## themselves, so they score as "flat", get large weights, and the baseline
## climbs straight through them.
##
## Tying the window to the data's own feature scale instead: the lag at which the
## autocorrelation of the first difference first crosses zero is a robust
## estimate of the half-width of the narrowest recurring feature.
.ff_auto_window <- function(y, min_hw = 4L, max_frac = 0.05) {
  n <- length(y)
  d <- diff(y)
  d <- d - mean(d)
  max_lag <- max(10L, min(500L, n %/% 10L))
  a <- stats::acf(d, lag.max = max_lag, plot = FALSE)$acf[-1]
  k <- which(a > 0)
  ## first difference of a smooth peak is positive-then-negative, so its acf
  ## dips negative and returns; the return crossing marks the feature scale
  first_pos_after_neg <- if (any(a < 0) && length(k)) k[k > which(a < 0)[1L]][1L] else NA_integer_
  hw <- if (is.finite(first_pos_after_neg)) first_pos_after_neg else max(min_hw, which.min(a))
  max(min_hw, min(as.integer(hw), as.integer(max_frac * n)))
}
# =============================================================================
# [NEW in v0.6.0] Whittaker-family asymmetric baselines: PSALSA, arPLS, AsLS
# =============================================================================
#
# All three minimise the same penalised least-squares criterion as FlatFit,
#
#     (y - b)' W (y - b)  +  lambda * b' D2' D2 b ,
#
# and differ only in how the weights W are updated from the residual
# d = y - b (peaks pointing up):
#
#   AsLS   (Eilers & Boelens 2005)  w = p      if d > 0,  1 - p  otherwise.
#   arPLS  (Baek et al. 2015)       w = 1 / (1 + exp(2 (d - (2 s - m)) / s)),
#                                   m, s = mean and sd of the negative d.
#   PSALSA (Oller-Moreno et al. 2014, "peaked signal's asymmetric least
#           squares"; the pybaselines psalsa() algorithm)
#                                   w = p * exp(-d / k)  if d > 0,  1 - p otherwise.
#          The exponential decay makes the weight of a point fall smoothly
#          with its height above the baseline, so tall peaks have practically
#          no pull while points within ~k of the baseline still count. This
#          is what keeps PSALSA from cutting into the foot of broad peaks.
#
# Automatic hyperparameters (all recorded in the result):
#   lambda  The second-difference Whittaker smoother has frequency response
#           1 / (1 + lambda (2 - 2 cos w)^2) ~ 1 / (1 + lambda w^4); its
#           half-power period is 2 pi lambda^(1/4) points. lambda is set so
#           that this period equals `lambda_scale` times the base width of
#           the widest peak (cluster) in the run,
#               lambda = (lambda_scale * W / (2 pi))^4 ,
#           with W = 2 * hw from MILE's peak-width estimator (estimate_hw():
#           widest contiguous run more than 3 noise sd above a coarse
#           morphological floor). In words: the baseline cannot bend on a
#           time scale shorter than lambda_scale peak widths. Because W is in
#           points, lambda is automatically invariant to the sampling rate.
#   k       (PSALSA) 3 x the noise sd (robust, from second differences):
#           points more than ~3 sd above the baseline are treated as peak.
#           (pybaselines' default k = sd(y)/10 is available as k = "std10".)
#   p       0.5 for PSALSA (pybaselines default), 0.01 for AsLS.
# Any of them can be given explicitly through `whittaker_args`.

## Robust noise sd from second differences (same estimator as MILE).
.ff_noise_sigma <- function(y) {
  d2 <- diff(y, differences = 2L)
  stats::mad(d2) / sqrt(6)
}

## Whittaker solve (W + lambda D2'D2) b = W y for any weights. The O(n)
## pentadiagonal Cholesky above is used directly: for chromatogram lengths it
## is faster than a general sparse solver (0.27 s vs 0.8 s at n = 24 001).
.ws_solve <- function(y, w, lambda) {
  n <- length(y)
  if (n < 5L) return(rep(stats::weighted.mean(y, w), n))
  .ff_whittaker(y, w, lambda)
}

#' Automatic Whittaker penalty from the widest peak (see block comment above).
whittaker_lambda_auto <- function(y, lambda_scale = 4, hw = NULL) {
  if (is.null(hw)) hw <- estimate_hw(y)
  W <- 2 * max(2L, as.integer(hw))
  (lambda_scale * W / (2 * pi))^4
}

#' Whittaker-family asymmetric baseline (PSALSA, arPLS or AsLS).
#'
#' @param y,x        signal and abscissa (x only checked for length).
#' @param method     "psalsa" (default), "arpls" or "asls".
#' @param lambda     smoothness penalty, sample units. NULL = automatic.
#' @param lambda_scale baseline cut-off period in widest-peak widths (auto
#'                   lambda only). Larger = stiffer baseline. Default 4.
#' @param hw         widest-peak half-width (points) for the auto lambda.
#'                   NULL = estimate_hw().
#' @param p          asymmetry. NULL = 0.5 (psalsa) / 0.01 (asls).
#' @param k          PSALSA decay height. NULL = 3 x noise sd; "std10" =
#'                   pybaselines default sd(y) / 10.
#' @param max_iter, tol  iteration limit and relative weight-change tolerance.
#' @param direction  "up" for peaks above the baseline, "down" for dips.
#' @param clamp      enforce baseline <= signal (>= for "down").
#' @return list(baseline, corrected, weights, lambda, lambda_scale, hw, p, k,
#'              iterations, converged, method)
whittaker_baseline <- function(
    y,
    x            = seq_along(y),
    method       = c("psalsa", "arpls", "asls"),
    lambda       = NULL,
    lambda_scale = 4,
    hw           = NULL,
    p            = NULL,
    k            = NULL,
    max_iter     = 50L,
    tol          = 1e-3,
    direction    = c("up", "down"),
    clamp        = FALSE
  ) {
  method    <- match.arg(method)
  direction <- match.arg(direction)
  stopifnot(length(x) == length(y), all(is.finite(y)))
  n <- length(y)
  if (n < 5L) stop("whittaker_baseline() needs at least 5 points.")

  sgn <- if (direction == "up") 1 else -1
  z   <- sgn * y

  if (is.null(hw)) hw <- estimate_hw(z)
  lambda_auto <- is.null(lambda)
  if (lambda_auto) lambda <- whittaker_lambda_auto(z, lambda_scale = lambda_scale, hw = hw)

  if (is.null(p)) p <- if (method == "asls") 0.01 else 0.5
  sig <- .ff_noise_sigma(z)
  if (!is.finite(sig) || sig <= 0) sig <- max(stats::sd(z), 1) * 1e-6
  if (method == "psalsa") {
    if (is.null(k)) {
      k <- 3 * sig
    } else if (identical(k, "std10")) {
      k <- stats::sd(z) / 10
    }
    k <- max(as.numeric(k), .Machine$double.eps)
  } else {
    k <- NA_real_
  }

  w    <- rep(1, n)
  b    <- .ws_solve(z, w, lambda)
  it   <- 0L
  conv <- FALSE
  for (i in seq_len(max_iter)) {
    d <- z - b
    w_new <- switch(method,
      asls   = ifelse(d > 0, p, 1 - p),
      psalsa = ifelse(d > 0, p * exp(-pmin(d / k, 700)), 1 - p),
      arpls  = {
        dn <- d[d < 0]
        if (length(dn) < 3L) { w } else {
          mm <- mean(dn); ss <- stats::sd(dn)
          if (!is.finite(ss) || ss <= 0) w else
            1 / (1 + exp(pmin(2 * (d - (2 * ss - mm)) / ss, 700)))
        }
      }
    )
    chg <- sqrt(sum((w_new - w)^2)) / max(sqrt(sum(w^2)), .Machine$double.eps)
    w   <- w_new
    b   <- .ws_solve(z, w, lambda)
    it  <- i
    if (chg < tol) { conv <- TRUE; break }
  }

  if (clamp) b <- pmin(b, z)
  b <- sgn * b
  list(
    baseline     = b,
    corrected    = y - b,
    weights      = w,
    lambda       = lambda,
    lambda_auto  = lambda_auto,
    lambda_scale = lambda_scale,
    hw           = as.integer(hw),
    p            = p,
    k            = k,
    iterations   = it,
    converged    = conv,
    method       = method
  )
}

psalsa_baseline <- function(y, x = seq_along(y), ...) whittaker_baseline(y, x, method = "psalsa", ...)
arpls_baseline  <- function(y, x = seq_along(y), ...) whittaker_baseline(y, x, method = "arpls",  ...)
asls_baseline   <- function(y, x = seq_along(y), ...) whittaker_baseline(y, x, method = "asls",   ...)

## ---------------------------------------------------------------------------
##  MILE: baseline correction using Morphological and Iterative Local Extremum
##
##  [MODIFIED v0.6.0] Documentation of the automatic structuring-element
##  width only (see estimate_hw()); no change in behaviour. estimate_hw() is
##  also used by the new PSALSA / arPLS / AsLS baselines to set lambda.
##
##  Reconstruction of the algorithm described in:
##
##    Ze-Yin Dong & Zi-Hang Yu (2023). "Baseline correction using morphological
##    and iterative local extremum (MILE)". Chemometrics and Intelligent
##    Laboratory Systems 240, 104908. doi:10.1016/j.chemolab.2023.104908
##
##  The algorithm follows four steps:
##
##    1. find all local extrema of the measured spectrum by derivation
##    2. estimate a coarse baseline by PCHIP interpolation of those extrema
##    3. find the local extrema of the coarse baseline again and update them
##       from their adjacent data; iterate the update several times
##    4. subtract the estimated baseline
##
##  Lineage: the "morphology + iterate + smooth" algorithm has its roots in
##  Perez-Pueyo et al. (2010) and Koch et al., J. Raman Spectrosc. 48:336-342
##  (2017) (IMM / "mormol").
## ---------------------------------------------------------------------------


## --- Helper functions --------------------------------------------------------

## Reflect-pad a vector by `hw` points on each side.
.pad_reflect <- function(y, hw) {
  n <- length(y)
  if (hw < 1) return(y)
  hw  <- min(hw, n - 1L)
  lft <- rev(y[2:(hw + 1L)])
  rgt <- rev(y[(n - hw):(n - 1L)])
  c(lft, y, rgt)
}

## Rolling min / max over a flat structuring element of half-width hw.
.roll <- function(y, hw, what = c("min", "max")) {
  what <- match.arg(what)
  if (hw < 1) return(y)
  n   <- length(y)
  hw  <- min(hw, n - 1L)
  yp  <- .pad_reflect(y, hw)
  out <- yp[seq_len(n)]
  f   <- if (what == "min") pmin else pmax
  for (k in seq_len(2L * hw)) out <- f(out, yp[(1L + k):(n + k)])
  out
}

erosion  <- function(y, hw) .roll(y, hw, "min")
dilation <- function(y, hw) .roll(y, hw, "max")
## Opening = dilation(erosion(.)): removes positive features narrower than the
## structuring element while leaving the underlying trend essentially in place.
opening  <- function(y, hw) dilation(erosion(y, hw), hw)
closing  <- function(y, hw) erosion(dilation(y, hw), hw)

## Moving average (reflect-padded), used as the "derivation" pre-smoother.
.movavg <- function(y, hw) {
  if (hw < 1) return(y)
  n  <- length(y)
  yp <- .pad_reflect(y, hw)
  cs <- cumsum(c(0, yp))
  w  <- 2L * min(hw, n - 1L) + 1L
  (cs[(w + 1L):(n + w)] - cs[seq_len(n)]) / w
}

## Mollifier (bump function) kernel, as in Koch et al. 2017.
.mollifier <- function(hw) {
  if (hw < 1) return(1)
  u <- seq(-1, 1, length.out = 2L * hw + 1L)
  k <- numeric(length(u))
  inside    <- abs(u) < 1
  k[inside] <- exp(-1 / (1 - u[inside]^2))
  k / sum(k)
}

.mollify <- function(y, hw) {
  if (hw < 1) return(y)
  k   <- .mollifier(hw)
  n   <- length(y)
  yp  <- .pad_reflect(y, hw)
  out <- numeric(n)
  for (j in seq_along(k)) out <- out + k[j] * yp[j:(j + n - 1L)]
  out
}

## Robust noise sigma from second differences.
noise_sigma <- function(y) {
  d2 <- diff(y, differences = 2L)
  stats::mad(d2) / sqrt(6)
}

## Indices of strict local minima, found from the sign of the first
## difference. Plateaus are handled by forward-filling zero signs.
local_minima <- function(y) {
  d  <- diff(y)
  s  <- sign(d)
  nz <- which(s != 0)
  if (length(nz) < 2L) return(integer(0))
  idx <- cumsum(s != 0)
  idx[idx == 0L] <- 1L
  sf <- s[nz][idx]                      # forward-filled signs
  which(sf[-length(sf)] < 0 & sf[-1] > 0) + 1L
}

local_maxima <- function(y) local_minima(-y)

## --- PCHIP -------------------------------------------------------------------

## Piecewise Cubic Hermite Interpolating Polynomial, Fritsch-Carlson slopes.
## Same construction as MATLAB's pchip() / SciPy's PchipInterpolator:
## shape-preserving, no overshoot between knots.
pchip <- function(xi, yi, xout) {
  o <- order(xi); xi <- xi[o]; yi <- yi[o]
  keep <- c(TRUE, diff(xi) > 0)
  xi <- xi[keep]; yi <- yi[keep]
  n <- length(xi)
  if (n == 1L) return(rep(yi, length(xout)))
  if (n == 2L) return(yi[1] + (yi[2] - yi[1]) * (xout - xi[1]) / (xi[2] - xi[1]))
  
  h  <- diff(xi)
  dl <- diff(yi) / h                    # secant slopes
  d  <- numeric(n)
  
  ## interior: weighted harmonic mean, zero at sign changes / flats
  if (n > 2L) {
    k  <- 2:(n - 1L)
    dm <- dl[k - 1L]; dp <- dl[k]
    hm <- h[k - 1L];  hp <- h[k]
    ok <- dm * dp > 0
    w1 <- 2 * hp + hm
    w2 <- hp + 2 * hm
    di <- numeric(length(k))
    di[ok] <- (w1[ok] + w2[ok]) / (w1[ok] / dm[ok] + w2[ok] / dp[ok])
    d[k] <- di
  }
  
  ## endpoints: one-sided three-point rule with the usual shape guards
  .end <- function(h1, h2, d1, d2) {
    de <- ((2 * h1 + h2) * d1 - h1 * d2) / (h1 + h2)
    if (sign(de) != sign(d1)) {
      de <- 0
    } else if (sign(d1) != sign(d2) && abs(de) > abs(3 * d1)) {
      de <- 3 * d1
    }
    de
  }
  if (n > 2L) {
    d[1] <- .end(h[1], h[2], dl[1], dl[2])
    d[n] <- .end(h[n - 1L], h[n - 2L], dl[n - 1L], dl[n - 2L])
  } else {
    d[1] <- dl[1]; d[n] <- dl[1]
  }
  
  ## Hermite evaluation
  j   <- findInterval(xout, xi, all.inside = TRUE)
  hj  <- h[j]
  t   <- (xout - xi[j]) / hj
  t2  <- t * t; t3 <- t2 * t
  h00 <-  2 * t3 - 3 * t2 + 1
  h10 <-      t3 - 2 * t2 + t
  h01 <- -2 * t3 + 3 * t2
  h11 <-      t3 -     t2
  h00 * yi[j] + h10 * hj * d[j] + h01 * yi[j + 1L] + h11 * hj * d[j + 1L]
}

## --- Automatic structuring-element size --------------------------------------

## Two ways to pick the structuring element automatically.
##
## [MODIFIED v0.6.0, documentation only] How the default is chosen, in full:
##   hw is the HALF-width (points) of the flat structuring element. MILE needs
##   hw >= half the base width of the widest genuine peak; any smaller and
##   the opening rides up into the peak. The default "peakwidth" estimator
##   (1) builds a deliberately oversized floor, opening(y, n/10) smoothed by a
##       moving average of the same half-width;
##   (2) measures the noise sd robustly from second differences
##       (noise_sigma(): MAD(diff2) / sqrt(6));
##   (3) finds every contiguous run of points more than k = 3 noise sd above
##       that floor -- each run is one peak or one cluster of merged peaks at
##       its base;
##   (4) returns hw = ceiling(quantile(run lengths, q) / 2) with q = 1 (the
##       widest run). Nothing is fitted; the only constants are n/10 (the
##       floor must be wider than any peak), k = 3 sd and q = 1. The value
##       used is stored per run as mile_hw. Pass mile_args = list(hw = ...) to
##       fix it.
##
## "peakwidth" (default): flatten the signal with a deliberately oversized
##   opening (half-width n/frac), so that *every* real peak shows up in the
##   residual, then measure the width of the contiguous runs that stand more
##   than k sigma above that floor and take half the widest. Assumes the
##   widest peak is narrower than n/frac and that the baseline has no
##   structure on that scale.
##
## "stability": grow hw until the opening stops changing. Cheaper to justify
##   but it systematically undershoots when broad and narrow features coexist.
estimate_hw <- function(y, method = c("peakwidth", "stability"), ...) {
  method <- match.arg(method)
  if (method == "peakwidth") .hw_peakwidth(y, ...) else .hw_stability(y, ...)
}

.hw_peakwidth <- function(y, frac = 10, k = 3, q = 1.0) {
  n   <- length(y)
  hw0 <- max(4L, floor(n / frac))
  b0  <- .movavg(opening(y, hw0), hw0)
  sg  <- noise_sigma(y)
  if (!is.finite(sg) || sg <= 0) sg <- .Machine$double.eps
  rr <- rle((y - b0) > k * sg)
  w  <- rr$lengths[rr$values]
  if (!length(w)) return(2L)
  max(2L, ceiling(as.numeric(stats::quantile(w, q)) / 2))
}

## Grow the structuring element until the opening stops changing: while hw is
## smaller than the peaks, each increment eats further into them and the
## opening moves a lot; once hw exceeds the widest peak the opening just
## tracks the baseline and successive openings barely differ. Same stopping
## rule as pybaselines' optimize_window() / Chen & Dai (2018).
.hw_stability <- function(y, max_hw = NULL, tol = 3e-4, max_hits = 3L) {
  n <- length(y)
  if (is.null(max_hw)) max_hw <- max(4L, floor(n / 4))
  rng <- diff(range(y))
  if (rng <= 0) return(2L)
  prev <- y
  hits <- 0L
  for (hw in seq_len(max_hw)) {
    o     <- opening(y, hw)
    delta <- mean(abs(o - prev)) / rng
    prev  <- o
    if (delta < tol) {
      hits <- hits + 1L
      if (hits >= max_hits) return(max(2L, hw - max_hits + 1L))
    } else hits <- 0L
  }
  max_hw
}

## --- MILE --------------------------------------------------------------------

#' Baseline correction using morphological and iterative local extremum
#'
#' @param y          numeric vector, the measured signal.
#' @param x          optional abscissa (defaults to seq_along(y)). Only used
#'                   for the PCHIP knot geometry; may be irregular.
#' @param hw         half-width of the flat structuring element, in points.
#'                   Should be at least the half-width of the widest genuine
#'                   peak. NULL = estimate from the data.
#' @param n_iter     maximum number of local-extremum update sweeps (step 3).
#' @param direction  "up" for emission-type peaks sitting on the baseline
#'                   (baseline = lower envelope, knots = local minima);
#'                   "down" for absorption-type dips (baseline = upper
#'                   envelope, knots = local maxima).
#' @param smooth_hw  half-width of the moving average applied before extremum
#'                   detection. NULL = hw %/% 4, which is usually enough to
#'                   stop noise from generating thousands of spurious knots.
#' @param k_morph    morphological screening threshold, in units of sigma.
#'                   A candidate knot is discarded if it sits more than
#'                   k_morph * sigma above the morphological opening -- i.e.
#'                   if it is a saddle between overlapping peaks rather than a
#'                   true return to baseline. Set to Inf to disable screening.
#' @param k_iter     update threshold, in units of sigma, for step 3. A knot is
#'                   only pulled down if it exceeds the chord through its two
#'                   neighbours by more than k_iter * sigma. This is what stops
#'                   the iteration from flattening genuinely curved baselines.
#' @param relax      relaxation factor in (0, 1] for each update sweep.
#' @param max_gap    a knot is only eligible for the iterative update if both
#'                   its neighbours are within this many points. NULL = Inf
#'                   (no restriction); set it to a few times hw if you want the
#'                   update confined to densely-knotted, overlapping regions.
#' @param clamp_to_opening  if TRUE, the iterative update may not push a knot
#'                   below the morphological opening. This protects genuine
#'                   baseline curvature from being flattened, at the cost of
#'                   making the result more sensitive to an undersized hw
#'                   (the opening itself rides up the peaks when hw is small).
#' @param mollify_hw half-width of an optional final mollifier smoothing of the
#'                   baseline. 0 = off.
#' @param clamp      if TRUE, enforce baseline <= signal (>= for "down").
#' @param knot_stat  how the ordinate of each knot is read off the data:
#'                   "smooth" (value of the pre-smoothed signal), "min" or
#'                   "median" of a +/- smooth_hw neighbourhood of the raw data.
#'
#' @return list with `baseline`, `corrected`, the knot set actually used,
#'         the parameters resolved, and the number of sweeps performed.
mile <- function(
    y,
    x                = seq_along(y),
    hw               = NULL,
    n_iter           = 20L,
    direction        = c("up", "down"),
    smooth_hw        = NULL,
    k_morph          = 3,
    k_iter           = 2,
    relax            = 1,
    max_gap          = NULL,
    clamp_to_opening = TRUE,
    mollify_hw       = 0L,
    clamp            = FALSE,
    knot_stat        = c("median", "smooth", "min")
) {
  
  direction <- match.arg(direction)
  stopifnot(length(x) == length(y), all(is.finite(y)))
  n <- length(y)
  
  ## Work internally with peaks pointing up; flip back at the end.
  sgn <- if (direction == "up") 1 else -1
  z   <- sgn * y
  
  sig <- noise_sigma(z)
  if (!is.finite(sig) || sig <= 0) sig <- .Machine$double.eps
  
  if (is.null(hw))        hw        <- estimate_hw(z)
  if (is.null(smooth_hw)) smooth_hw <- max(0L, hw %/% 4L)
  
  ## --- step 1: local extrema of the (lightly smoothed) spectrum ------------
  zs <- .movavg(z, smooth_hw)
  kn <- local_minima(zs)
  
  ## --- morphological screening ---------------------------------------------
  ## opening(zs, hw) is a lower envelope that has had every feature narrower
  ## than the structuring element removed. A local minimum lying well above it
  ## is a valley *between* overlapping peaks, not a point on the baseline.
  op <- opening(zs, hw)
  if (length(kn) && is.finite(k_morph)) {
    kn <- kn[(zs[kn] - op[kn]) <= k_morph * sig]
  }
  
  ## always anchor both ends
  kn <- sort(unique(c(1L, kn, n)))
  
  ## Knot ordinates. "smooth" reads the pre-smoothed signal (low noise bias,
  ## slight upward bias from peak flanks); "min" takes the raw minimum in a
  ## +/- smooth_hw neighbourhood (no flank bias, but biased down by roughly
  ## E[min of 2*smooth_hw+1 normals] * sigma).
  knot_stat <- match.arg(knot_stat)
  w  <- max(1L, smooth_hw)
  kv <- switch(knot_stat,
               smooth = zs[kn],
               min    = vapply(kn, function(i) min(z[max(1L, i - w):min(n, i + w)]), numeric(1)),
               median = vapply(kn, function(i) stats::median(z[max(1L, i - w):min(n, i + w)]), numeric(1)))
  
  ## --- step 2: coarse baseline by PCHIP -----------------------------------
  b <- pchip(x[kn], kv, x)
  
  ## --- step 3: iterative local-extremum update ----------------------------
  ## Local maxima of the coarse baseline can only arise where a knot was
  ## dragged upward by an unresolved peak. Pull such knots down toward the
  ## chord joining their neighbours, but only when the excess is well above
  ## the noise -- a smoothly curved baseline also produces knot-level maxima,
  ## and those must be preserved.
  max_gap_pts <- if (is.null(max_gap)) Inf else max_gap
  it <- 0L
  if (n_iter > 0L && length(kn) >= 3L) {
    for (sweep in seq_len(n_iter)) {
      m <- length(kv)
      i <- 2:(m - 1L)
      ## chord value at knot i from its two neighbours
      tt     <- (x[kn[i]] - x[kn[i - 1L]]) / (x[kn[i + 1L]] - x[kn[i - 1L]])
      chord  <- kv[i - 1L] + tt * (kv[i + 1L] - kv[i - 1L])
      excess <- kv[i] - chord
      is_max <- kv[i] > kv[i - 1L] & kv[i] > kv[i + 1L]
      ## Only knots in *dense* stretches are candidates. A saddle between
      ## overlapping peaks sits a peak-width away from its neighbours; a knot
      ## that is isolated is far more likely to be a real point on a curved
      ## baseline, and pulling it to the chord would flatten genuine curvature.
      gap <- pmax(kn[i] - kn[i - 1L], kn[i + 1L] - kn[i])
      hit <- which(is_max & excess > k_iter * sig & gap <= max_gap_pts)
      if (!length(hit)) break
      ## Never pull a knot below the morphological opening: the opening is a
      ## valid lower envelope of the data, so it bounds how far a knot can
      ## legitimately be wrong. This is what stops the iteration from eating
      ## real baseline curvature when max_gap is loose.
      upd <- kv[i[hit]] - relax * excess[hit]
      new <- if (clamp_to_opening) pmax(upd, op[kn[i[hit]]]) else upd
      if (max(abs(new - kv[i[hit]])) < 1e-12) break   # all candidates already
      kv[i[hit]] <- new                               # sit on the opening
      it <- sweep
    }
    b <- pchip(x[kn], kv, x)
  }
  
  ## --- optional smoothing + clamp -----------------------------------------
  if (mollify_hw > 0L) b <- .mollify(b, mollify_hw)
  if (clamp)           b <- pmin(b, z)
  
  ## --- step 4: subtract ----------------------------------------------------
  b <- sgn * b
  list(
    baseline   = b,
    corrected  = y - b,
    knots      = kn,
    knot_x     = x[kn],
    knot_y     = sgn * kv,
    sigma      = sig,
    hw         = hw,
    smooth_hw  = smooth_hw,
    iterations = it
  )
}
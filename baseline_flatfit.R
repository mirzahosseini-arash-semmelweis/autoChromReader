## ---------------------------------------------------------------------------
##  flatfit.R -- Flatness-weighted baseline fit.
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
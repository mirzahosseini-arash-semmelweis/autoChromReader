## ---------------------------------------------------------------------------
##  MILE: baseline correction using Morphological and Iterative Local Extremum
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
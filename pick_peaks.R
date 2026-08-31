## ---------------------------------------------------------------------------
##  Robust peak picking for (baseline-corrected) chromatograms
##
##  Three independent detectors, deliberately, because no single statistic
##  finds every kind of peak:
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
##  exponential-Gaussian-hybrid (EGH) profiles so that the apex of a genuinely
##  merged component can be recovered rather than guessed.
##
##  References
##    Du P., Kibbe W.A., Lin S.M. (2006) Bioinformatics 22:2059-2065.
##    Lan K., Jorgenson J.W. (2001) J. Chromatogr. A 915:1-13.  (EGH)
##    Tautenhahn R. et al. (2008) BMC Bioinformatics 9:504.     (centWave)
## ---------------------------------------------------------------------------


## --- Helper functions --------------------------------------------------------

## convolve y with a kernel of odd length, reflect-padded, O(n * length(k))
.convolve <- function(y, k) {
  hw  <- (length(k) - 1L) %/% 2L
  if (hw < 1L) return(y * sum(k))
  n   <- length(y)
  yp  <- .pad_reflect(y, hw)
  out <- numeric(n)
  for (j in seq_along(k)) out <- out + k[j] * yp[j:(j + n - 1L)]
  out
}

## Savitzky-Golay coefficients: polynomial order `ord`, derivative `d`.
.sg_coef <- function(hw, ord = 3L, d = 0L) {
  z <- -hw:hw
  X <- outer(z, 0:ord, "^")
  C <- solve(crossprod(X), t(X))
  factorial(d) * C[d + 1L, ]
}

sg_filter <- function(y, hw, ord = 3L, d = 0L, dx = 1) {
  if (hw < ord) hw <- ord + 1L
  .convolve(y, .sg_coef(hw, ord, d)) / dx^d
}

#' Noise measured the way a chromatographer measures it: in peak-free windows.
#' Splits the trace into windows, removes a linear trend from each, and takes a
#' low quantile of the per-window sd -- windows containing peaks have a large
#' sd and are excluded automatically. Also returns the pooled residual of the
#' quiet windows, which is the correct pool to bootstrap from (bootstrapping
#' the whole residual smuggles peak curvature into the null and inflates every
#' threshold derived from it).
noise_profile <- function(y, win = 201L, q = 0.25) {
  n     <- length(y)
  win   <- min(win, max(11L, n %/% 8L))
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
  pool <- unlist(resid[quiet], use.names = FALSE)
  list(
    sd      = sigma,
    pool    = pool,
    win_sd  = sds,
    n_quiet = length(quiet)
  )
}

## sd of a convolution coefficient under the empirical noise autocovariance --
## exact for any stationary noise, no simulation needed.
.coef_sd <- function(k, gam) {
  m <- length(k); L <- min(m - 1L, length(gam) - 1L)
  s <- gam[1] * sum(k^2)
  if (L >= 1L) for (lag in seq_len(L))
    s <- s + 2 * gam[lag + 1L] * sum(k[1:(m - lag)] * k[(1 + lag):m])
  sqrt(max(s, 0))
}

## Ricker / Mexican-hat wavelet, L2-normalised, at scale s (in points).
.ricker <- function(s) {
  hw <- max(1L, ceiling(4 * s))
  t  <- -hw:hw
  k  <- (2 / (sqrt(3 * s) * pi^0.25)) * (1 - (t / s)^2) * exp(-t^2 / (2 * s^2))
  k / sqrt(sum(k^2))
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
.surrogate <- function(r, block) {
  n   <- length(r)
  nb  <- ceiling(n / block)
  st  <- sample.int(n, nb, replace = TRUE)
  idx <- as.vector(vapply(st, function(s) ((s + 0:(block - 1L) - 1L) %% n) + 1L,
                          integer(block)))
  r[idx[seq_len(n)]]
}

## Autocorrelation length of the residual: first lag where acf drops below 1/e.
.acf_len <- function(r, max_lag = 50L) {
  a <- stats::acf(r, lag.max = max_lag, plot = FALSE)$acf[-1]
  w <- which(a < exp(-1))
  if (!length(w)) max_lag else w[1]
}

#' Largest prominence attributable to noise alone.
#' @param fwer probability that a noise-only chromatogram yields ANY peak above
#'   the returned cutoff. 0.05 => at most a 1-in-20 chance of a single false
#'   positive anywhere in the run.
calibrate_prominence <- function(y, nsim = 25L, fwer = 0.05, np = NULL) {
  if (is.null(np)) np <- noise_profile(y)
  r  <- np$pool
  b  <- max(2L, .acf_len(r) * 3L)
  n  <- length(y)
  mx <- vapply(seq_len(nsim), function(i) {
    s <- .surrogate(r, b)[seq_len(min(n, length(r)))]
    if (length(s) < n) s <- rep_len(s, n)
    max(prominence(s)$prominence)
  }, numeric(1))
  list(
    cutoff   = as.numeric(stats::quantile(mx, 1 - fwer)),
    sd       = np$sd,
    block    = b,
    null_max = mx,
    gamma = stats::acf(r, lag.max = 200L, type = "covariance",
                       plot = FALSE)$acf[, 1, 1]
  )
}

## --- CWT ridge detection -----------------------------------------------------

#' Ricker CWT plus ridge linking (Du et al. 2006).
#' @param scales wavelet scales in POINTS. Should bracket sigma of the
#'   narrowest and widest expected peak.
#' @return list(W = n x nscale coefficient matrix, ridges = data.frame)
cwt_ridges <- function(y, scales, gap_max = 3L, win_fac = 2) {
  n  <- length(y)
  W  <- vapply(scales, function(s) .convolve(y, .ricker(s)), numeric(n))
  ns <- length(scales)

  ## local maxima of the coefficient at each scale
  lmax <- lapply(seq_len(ns), function(j) {
    w <- W[, j]
    which(c(FALSE, diff(w) > 0)[seq_len(n)] & c(diff(w) < 0, FALSE))
  })

  ## link from coarse to fine, each ridge tracking its position downwards
  ridges <- list()
  active <- list()
  for (j in rev(seq_len(ns))) {
    cand <- lmax[[j]]
    used <- logical(length(cand))
    win  <- max(1L, ceiling(win_fac * scales[j]))
    if (length(active)) {
      drop <- logical(length(active))
      for (a in seq_along(active)) {
        p <- active[[a]]$pos
        linked <- FALSE
        if (length(cand)) {
          d <- abs(cand - p)
          ok <- which(!used & d <= win)
          if (length(ok)) {
            k <- ok[which.min(d[ok])]
            used[k] <- TRUE
            active[[a]]$pos <- cand[k]
            active[[a]]$idx <- c(active[[a]]$idx, cand[k])
            active[[a]]$sc  <- c(active[[a]]$sc, j)
            active[[a]]$gap <- 0L
            linked <- TRUE
          }
        }
        if (!linked) {
          active[[a]]$gap <- active[[a]]$gap + 1L
          if (active[[a]]$gap > gap_max) {
            ridges[[length(ridges) + 1L]] <- active[[a]]
            drop[a] <- TRUE
          }
        }
      }
      if (any(drop)) active <- active[!drop]
    }
    for (k in which(!used))
      active[[length(active) + 1L]] <- list(pos = cand[k], idx = cand[k], sc = j, gap = 0L)
  }
  ridges <- c(ridges, active)

  if (!length(ridges)) return(list(W = W, ridges = NULL, scales = scales))
  rid <- do.call(rbind, lapply(ridges, function(rg) {
    v <- W[cbind(rg$idx, scales[rg$sc] * 0 + rg$sc)]
    best <- which.max(v)
    data.frame(
      idx      = rg$idx[best],
      idx_fine = rg$idx[length(rg$idx)],
      scale    = scales[rg$sc[best]],
      len      = length(rg$sc),
      coef     = v[best]
    )
  }))
  list(
    W      = W,
    ridges = rid[order(rid$idx), ],
    scales = scales
  )
}

## --- EGH peak shape ----------------------------------------------------------

#' Exponential-Gaussian hybrid (Lan & Jorgenson 2001). Algebraic, so it is
#' stable to fit, and it models the tailing that every real HPLC peak has.
#' tau > 0 tails to the right, tau < 0 fronts.
egh <- function(t, H, tR, sg, tau) {
  den <- 2 * sg^2 + tau * (t - tR)
  out <- numeric(length(t))
  ok  <- den > 0
  out[ok] <- H * exp(-(t[ok] - tR)^2 / den[ok])
  out
}

egh_width <- function(sg, tau, frac = 0.5) {
  L <- log(1 / frac)
  d <- sqrt(L^2 * tau^2 + 8 * sg^2 * L)
  cbind(left = (L * tau - d) / 2, right = (L * tau + d) / 2, width = d)
}

egh_area <- function(H, tR, sg, tau, rel_tol = 1e-8) {
  span <- 25 * sg + 25 * abs(tau)
  a <- try(stats::integrate(function(t) egh(t, H, tR, sg, tau),
                            tR - span, tR + span,
                            rel.tol = rel_tol, subdivisions = 500L)$value,
           silent = TRUE)
  if (inherits(a, "try-error")) {            # fall back to a fine grid
    t <- seq(tR - span, tR + span, length.out = 20001L)
    a <- sum(egh(t, H, tR, sg, tau)) * (t[2] - t[1])
  }
  a
}

#' Fit a sum of EGH components to one cluster of merged peaks.
#' @param starts data.frame(tR, H, sg) initial guesses, one row per component.
fit_cluster <- function(
    x,
    y,
    starts,
    tau0      = 0.35,
    maxit     = 800L,
    n_restart = 4L
  ) {
  k   <- nrow(starts)
  tR0 <- starts$tR
  sg0 <- pmax(starts$sg, .Machine$double.eps)
  unpack <- function(th) list(tR  = tR0 + th[1:k] * sg0,
                              sg  = sg0 * exp(pmin(pmax(th[(k + 1):(2 * k)], -6), 6)),
                              tau = th[(2 * k + 1):(3 * k)] * sg0)
  basis <- function(q) vapply(seq_len(k), function(i)
    egh(x, 1, q$tR[i], q$sg[i], q$tau[i]), numeric(length(x)))
  
  solveH <- function(Phi) {
    H <- try(as.numeric(qr.solve(Phi, y)), silent = TRUE)
    if (inherits(H, "try-error") || any(!is.finite(H))) return(NULL)
    if (any(H < 0)) {
      keep <- H >= 0
      if (!any(keep)) return(rep(0, k))
      H2 <- rep(0, k)
      h  <- try(as.numeric(qr.solve(Phi[, keep, drop = FALSE], y)), silent = TRUE)
      if (inherits(h, "try-error")) return(NULL)
      H2[keep] <- pmax(h, 0)
      return(H2)
    }
    H
  }
  
  obj <- function(th) {
    q <- unpack(th)
    if (any(q$sg <= 0)) return(1e12)
    Phi <- basis(q)
    if (any(colSums(Phi) <= 0)) return(1e12)
    H <- solveH(Phi)
    if (is.null(H)) return(1e12)
    sum((Phi %*% H - y)^2)
  }
  
  th0  <- c(rep(0, k), rep(0, k), rep(tau0, k))
  best <- NULL
  for (r in seq_len(max(1L, n_restart))) {
    th <- if (r == 1L) th0 else th0 + stats::rnorm(length(th0), 0, 0.25)
    f <- stats::optim(th, obj, method = "Nelder-Mead",
                      control = list(maxit = maxit, reltol = 1e-12))
    f <- stats::optim(f$par, obj, method = "BFGS", control = list(maxit = maxit))
    if (is.null(best) || f$value < best$value) best <- f
  }
  fit <- best
  
  q    <- unpack(fit$par)
  q$H  <- solveH(basis(q))
  comp <- vapply(seq_len(k), function(i) egh(x, q$H[i], q$tR[i], q$sg[i], q$tau[i]),
                 numeric(length(x)))
  dxm  <- mean(diff(x))
  
  w50 <- unname(egh_width(q$sg, q$tau, 0.50))
  colnames(w50) <- c("left", "right", "width")
  
  data.frame(
    component   = seq_len(k),
    rt          = q$tR,
    tR          = q$tR,
    H           = q$H,
    height      = q$H,
    fwhm        = w50[, "width"],
    fwhm_left   = q$tR + w50[, "left"],
    fwhm_right  = q$tR + w50[, "right"],
    sigma       = q$sg,
    tau         = q$tau,
    area        = vapply(seq_len(k), function(i)
      egh_area(q$H[i], q$tR[i], q$sg[i], q$tau[i]), numeric(1)),
    area_window = colSums(comp) * dxm,
    rss         = fit$value,
    npar        = 4L * k,
    nobs        = length(x),
    converged   = fit$convergence == 0,
    row.names   = NULL
  )
}

## --- Main peak finder --------------------------------------------------------
#' Robust peak picking for a baseline-corrected chromatogram
#'
#' @param x,y     retention time and (baseline-corrected) signal.
#' @param fwer    family-wise false-positive rate for the prominence cutoff.
#' @param min_prom  optional hard override of the calibrated cutoff.
#' @param snr_min minimum apex height / noise sd. Applied on top of prominence.
#' @param scales  CWT scales in points. NULL = geometric grid inferred from the
#'                widths of the prominence-detected peaks.
#' @param min_ridge_len  a CWT ridge must persist over this many scales.
#' @param shoulder_snr   a shoulder-only candidate must reach this coefficient
#'                       S/N to be reported.
#' @param nsim    bootstrap replicates for the noise calibration.
#'
#' @return data.frame, one row per peak:
#'   rt, idx, height, prominence, snr, fwhm, left/right (integration bounds),
#'   area, type ("peak" | "shoulder"), cluster, resolution (Rs to the previous
#'   peak in the same cluster).
pick_peaks <- function(
    x,
    y,
    fwer          = 0.05,
    min_prom      = NULL,
    snr_min       = 20,
    scales        = NULL,
    min_ridge_len = 4L,
    shoulder_snr  = NULL,
    nsim          = 25L,
    nsim_cwt      = 5L,
    do_shoulders  = TRUE
  ) {

  n    <- length(y)
  dx   <- mean(diff(x))
  np   <- noise_profile(y)
  sd_n <- np$sd

  ## ---- (1)+(2) prominence, with a bootstrapped cutoff ---------------------
  cal  <- calibrate_prominence(y, nsim = nsim, fwer = fwer, np = np)
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
  if (is.null(scales)) {
    lo <- max(1, min(sig_pts) / 2); hi <- max(sig_pts) * 3
    scales <- unique(round(exp(seq(log(lo), log(max(hi, lo * 4)), length.out = 12)), 2))
  }

  ## ---- (3) CWT ridges ------------------------------------------------------
  sh <- NULL
  if (do_shoulders) {
    cw <- cwt_ridges(y, scales)
    if (!is.null(cw$ridges)) {
      rid <- cw$ridges[cw$ridges$len >= min_ridge_len, , drop = FALSE]
      ## Noise level of each coefficient, computed exactly from the empirical
      ## noise autocovariance rather than assuming white noise.
      csd <- vapply(rid$scale, function(s) .coef_sd(.ricker(s), cal$gamma),
                    numeric(1))
      rid$snr <- rid$coef / csd
      ## Threshold calibrated the same way as the prominence cutoff: run the
      ## noise surrogates through the identical CWT + ridge pipeline and take
      ## the largest ridge S/N noise produces. Any bias in .coef_sd cancels,
      ## and it accounts for the multiple testing across ~n/(2s) positions and
      ## every scale.
      if (is.null(shoulder_snr)) {
        b  <- cal$block
        mx <- vapply(seq_len(nsim_cwt), function(i) {
          sg <- rep_len(.surrogate(np$pool, b), n)
          cs <- cwt_ridges(sg, scales)
          if (is.null(cs$ridges)) return(0)
          rr <- cs$ridges[cs$ridges$len >= min_ridge_len, , drop = FALSE]
          if (!nrow(rr)) return(0)
          max(rr$coef / vapply(rr$scale,
                function(s) .coef_sd(.ricker(s), cal$gamma), numeric(1)))
        }, numeric(1))
        shoulder_snr <- max(mx)
      }
      rid <- rid[rid$snr >= shoulder_snr, , drop = FALSE]
      ## drop ridges that merely re-detect an already-found apex
      if (nrow(rid) && nrow(pk)) {
        tol  <- pmax(3, 1.5 * rid$scale)
        dmin <- vapply(seq_len(nrow(rid)),
                       function(j) min(abs(rid$idx_fine[j] - pk$idx)), numeric(1))
        rid  <- rid[dmin > tol, , drop = FALSE]
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
  if (!nrow(out)) return(out)
  out <- out[order(out$idx), ]
  ## edge guard: prominence is ill-defined against a truncated flank
  edge <- max(3L, ceiling(min(scales)))
  out  <- out[out$idx > edge & out$idx <= n - edge, , drop = FALSE]
  if (!nrow(out)) return(out)

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

  ## Flat-topped (saturated) apexes: a clipped peak has a plateau, which puts
  ## a curvature maximum at each plateau EDGE and would otherwise be reported
  ## as two shoulders. Flag them and drop shoulders that fall inside.
  out$saturated <- vapply(seq_len(m), function(j) {
    i  <- out$idx[j]
    fl <- sum(y[out$left[j]:out$right[j]] >= 0.995 * y[i])
    fl > max(3, 0.25 * out$fwhm[j] / dx)
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
  attr(out, "noise_sd")     <- sd_n
  attr(out, "prom_cutoff")  <- cal$cutoff
  attr(out, "scales")       <- scales
  attr(out, "null_max")     <- cal$null_max
  attr(out, "shoulder_snr") <- if (do_shoulders) shoulder_snr else NA
  out
}

#' How many components does this cluster actually need?
#'
#' Fits k = 1..k_max EGH profiles and selects k by BIC. Extra components are
#' seeded at the largest positive residual of the previous fit, so a merged
#' peak that is invisible in both the trace AND its second derivative can still
#' be recovered -- which matters, because a shoulder at Rs ~ 0.45 with a 3:1
#' height ratio produces no separate maximum in y and none in y'' either.
#'
#' The BIC uses an EFFECTIVE sample size m / L, where L is the noise
#' autocorrelation length. Chromatographic noise is not white; using the raw
#' point count would make BIC think it has far more evidence than it does and
#' it would happily fit components to wiggles.
select_components <- function(
    x,
    y,
    seeds,
    k_max         = 3L,
    acl           = 1,
    bic_margin    = 10,
    min_area_frac = 0.01
  ) {
  m     <- length(y)
  m_eff <- max(5, m / max(1, acl))
  best  <- NULL; best_bic <- Inf; hist <- NULL
  st    <- seeds[order(-seeds$H), , drop = FALSE]
  for (k in seq_len(k_max)) {
    if (k <= nrow(st)) {
      s0 <- st[seq_len(k), , drop = FALSE]
    } else {
      if (is.null(best)) break
      prev <- best
      pred <- numeric(m)
      for (i in seq_len(nrow(prev)))
        pred <- pred + egh(x, prev$H[i], prev$tR[i], prev$sigma[i], prev$tau[i])
      res <- y - pred
      j <- which.max(res)
      if (res[j] <= 0) break
      s0 <- rbind(data.frame(tR = prev$tR, H = prev$height, sg = prev$sigma),
                  data.frame(tR = x[j],    H = res[j],      sg = stats::median(prev$sigma)))
    }
    f <- try(fit_cluster(x, y, s0), silent = TRUE)
    if (inherits(f, "try-error")) next
    ## a component carrying a negligible share of the cluster is a fitting
    ## artefact, not an analyte -- reject the whole model rather than report it
    if (k > 1L && min(f$area) / sum(f$area) < min_area_frac) next
    bic  <- m_eff * log(f$rss[1] / m) + 4 * k * log(m_eff)
    hist <- rbind(hist, data.frame(k = k, rss = f$rss[1], bic = bic))
    ## accept a richer model only on strong evidence (Kass & Raftery: a BIC
    ## drop of >10 is "very strong"). Without this the fit happily adds a
    ## component to every bit of tailing it cannot otherwise explain.
    if (bic < best_bic - (if (is.null(best)) 0 else bic_margin)) {
      best_bic <- bic; best <- f
    }
  }
  list(fit = best, bic = hist)
}

#' Refit every multi-component cluster with a sum of EGH profiles.
#' Returns the same peak table with apex, height and area replaced by the
#' fitted values where a fit succeeded.
deconvolve_clusters <- function(
    x,
    y,
    pk,
    min_members = 1L,
    pad         = 1.5,
    k_max       = 3L,
    min_snr     = 20,
    acl         = NULL
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
                                 k_max = max(k_max, nrow(st) + 1L), acl = acl),
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
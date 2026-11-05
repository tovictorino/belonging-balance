## psfm_gtre_zh() with a choice of Halton draw scheme:
##   draws = "shared"  : one set of draws for all individuals (as in psfm_gtre_zh.R, sfa 1.0.4)
##   draws = "perfirm" : each individual gets its own block of draws (as in sfa 1.2.0)
## Sourcing this file replaces psfm_gtre_zh(); with draws = "shared" it
## reproduces psfm_gtre_zh.R exactly. Licensed GPL (>= 2).

source("psfm_gtre_zh.R")

.make_draws_perfirm <- function(N, R, discard = 1000L, seed = NULL, clamp = 1e-6) {
  if (!requireNamespace("randtoolbox", quietly = TRUE))
    stop("draws = \"perfirm\" needs the 'randtoolbox' package")
  U <- randtoolbox::halton(N * R + discard, dim = 2, start = 1, normal = FALSE)
  U <- U[-seq_len(discard), , drop = FALSE]
  if (!is.null(seed)) {
    set.seed(seed)
    sh <- stats::runif(2L)
    U[, 1] <- (U[, 1] + sh[1]) %% 1
    U[, 2] <- (U[, 2] + sh[2]) %% 1
  }
  U[] <- pmin(pmax(U, clamp), 1 - clamp)
  list(r = matrix(stats::qnorm(U[, 1]), N, R, byrow = TRUE),
       h = matrix(stats::qnorm((U[, 2] + 1) / 2), N, R, byrow = TRUE))
}

## zr, zh: length-R vectors (shared draws) or N x R matrices (per-firm draws).
.zh_make_objective <- function(Yv, X, Zu, W, gid, zr, zh, s) {
  pfd <- is.matrix(zr)
  n  <- length(Yv); N <- nrow(W); R <- if (pfd) ncol(zr) else length(zr)
  kx <- ncol(X); pu <- ncol(Zu); qh <- ncol(W)
  npar <- 2L + kx + pu + qh
  i_b <- 3L:(2L + kx)
  i_d <- (3L + kx):(2L + kx + pu)
  i_g <- (3L + kx + pu):(2L + kx + pu + qh)
  CLIP <- .zh_constants$LIN_CLIP
  BIGN <- .zh_constants$BIG / N
  logR <- log(R)
  LOG_SQRT_2PI <- 0.5 * log(2 * pi)
  if (pfd) {
    zr_n <- zr[gid, , drop = FALSE]
    zh_n <- zh[gid, , drop = FALSE]
  } else {
    zr_n <- rep(zr, each = n)
  }

  function(p, deriv = FALSE) {
    sigv <- p[1]; sigr <- p[2]
    if (!all(is.finite(p)) || sigv <= 0 || sigr <= 0)
      return(list(value = .zh_constants$BIG, grad = rep(0, npar)))
    beta <- p[i_b]; delta <- p[i_d]; gam <- p[i_g]

    lh_raw <- as.numeric(W %*% gam)
    sh     <- sqrt(exp(pmin(pmax(lh_raw, -CLIP), CLIP)))
    lu_raw <- as.numeric(Zu %*% delta)
    su2    <- exp(pmin(pmax(lu_raw, -CLIP), CLIP))
    sg     <- sqrt(sigv^2 + su2)
    lam    <- sqrt(su2) / sigv
    base   <- s * (Yv - as.numeric(X %*% beta))
    E   <- (if (pfd) sh[gid] * zh_n else outer(sh[gid], zh)) + base - s * sigr * zr_n
    A   <- E / sg
    B   <- -A * lam
    lPB <- stats::pnorm(B, log.p = TRUE)
    lp  <- (log(2) - LOG_SQRT_2PI) - log(sg) - 0.5 * A * A + lPB
    S   <- rowsum(lp, gid)
    mx  <- S[cbind(seq_len(N), max.col(S, ties.method = "first"))]
    ES  <- exp(S - mx)
    sE  <- rowSums(ES)
    v   <- -(mx + log(sE) - logR)
    bad <- !is.finite(v)
    v[bad] <- BIGN
    value <- sum(v)
    if (!deriv) return(list(value = value))

    Wn  <- (ES / sE)[gid, , drop = FALSE]
    M   <- exp(-0.5 * B * B - LOG_SQRT_2PI - lPB)
    GE  <- Wn * (-(A + M * lam) / sg)
    gSG <- rowSums(Wn * (A^2 - 1 - M * B)) / sg
    gLM <- -rowSums(Wn * M * A)
    if (any(bad)) {
      ob <- bad[gid]
      GE[ob, ] <- 0; gSG[ob] <- 0; gLM[ob] <- 0
    }
    GEh  <- if (pfd) rowSums(GE * zh_n) else as.numeric(GE %*% zh)
    d_lh <- -(sh / 2) * as.numeric(rowsum(GEh, gid)) * (abs(lh_raw) < CLIP)
    d_lu <- -(gSG * su2 / (2 * sg) + gLM * lam / 2) * (abs(lu_raw) < CLIP)
    grad <- c(-sum(gSG * sigv / sg - gLM * lam / sigv),
              s * (if (pfd) sum(GE * zr_n) else sum(GE %*% zr)),
              s * as.numeric(crossprod(X, rowSums(GE))),
              as.numeric(crossprod(Zu, d_lu)),
              as.numeric(crossprod(W, d_lh)))
    list(value = value, grad = grad)
  }
}

psfm_gtre_zh <- function(formula, data, individual,
                         inefdec        = TRUE,
                         halton_num     = NULL,
                         halton_discard = 100L,
                         draws          = c("shared", "perfirm"),
                         intercept_zu   = TRUE,
                         intercept_zh   = TRUE,
                         start_val      = NULL,
                         maxit.bobyqa   = 300,
                         maxit.psoptim  = 10,
                         maxit.optim    = 200,
                         PSopt          = FALSE,
                         optHessian     = TRUE,
                         Method         = "L-BFGS-B",
                         analytic_grad  = TRUE,
                         compute_te     = TRUE,
                         rand.gtre      = NULL,
                         rand.psoptim   = NULL,
                         verbose        = FALSE) {

  call  <- match.call()
  t0    <- Sys.time()
  s     <- if (isTRUE(inefdec)) 1 else -1
  draws <- match.arg(draws)

  pf <- .parse_formula3(formula)
  if (!pf$y %in% names(data))       stop("response '", pf$y, "' not in data")
  if (!individual %in% names(data)) stop("individual '", individual, "' not in data")

  data <- as.data.frame(data)
  used <- unique(c(pf$y, individual, all.vars(pf$fx),
                   all.vars(pf$fzu), all.vars(pf$fzh)))
  cc   <- stats::complete.cases(data[, used, drop = FALSE])
  if (any(!cc)) {
    message(sum(!cc), " incomplete rows dropped.")
    data <- data[cc, , drop = FALSE]
  }
  data <- data[order(data[[individual]]), , drop = FALSE]

  Yv <- as.numeric(data[[pf$y]])
  X  <- stats::model.matrix(pf$fx,  data)
  Zu <- stats::model.matrix(pf$fzu, data)
  Zh <- stats::model.matrix(pf$fzh, data)
  if (!intercept_zu) Zu <- Zu[, colnames(Zu) != "(Intercept)", drop = FALSE]
  if (!intercept_zh) Zh <- Zh[, colnames(Zh) != "(Intercept)", drop = FALSE]
  if (ncol(Zu) == 0L) stop("transient margin has no columns; use '| 1 |' or intercept_zu = TRUE")
  if (ncol(Zh) == 0L) stop("persistent margin has no columns; use '| 1' or intercept_zh = TRUE")

  id    <- data[[individual]]
  indiv <- unique(id)
  N     <- length(indiv)
  idx   <- split(seq_len(nrow(data)), factor(id, levels = indiv))
  t_i   <- vapply(idx, length, integer(1))
  n     <- nrow(data)
  kx <- ncol(X); pu <- ncol(Zu); qh <- ncol(Zh)

  W <- do.call(rbind, lapply(idx, function(r) colMeans(Zh[r, , drop = FALSE])))
  colnames(W) <- colnames(Zh)
  nonint <- setdiff(colnames(Zh), "(Intercept)")
  if (length(nonint)) {
    wv <- vapply(nonint, function(j)
      max(vapply(idx, function(r)
        if (length(r) > 1) stats::sd(Zh[r, j]) else 0, numeric(1))), numeric(1))
    if (any(wv > 1e-8))
      warning("persistent-margin variable(s) vary within individuals: ",
              paste(names(wv)[wv > 1e-8], collapse = ", "),
              ". Individual means were used.")
  }

  npar   <- 2L + kx + pu + qh
  pnames <- c("sigv", "sigr", colnames(X),
              paste0("Zu:", colnames(Zu)), paste0("Zh:", colnames(W)))

  if (is.numeric(start_val)) {
    if (length(start_val) != npar)
      stop("start_val must have length ", npar, " (order: sigv, sigr, beta, delta, gamma)")
    start_v <- start_val
  } else {
    ols <- stats::lm.fit(X, Yv)
    b0  <- ols$coefficients; b0[is.na(b0)] <- 0
    res <- Yv - as.numeric(X %*% b0)
    am  <- vapply(idx, function(r) mean(res[r]), numeric(1))
    ew  <- res - rep(am, times = t_i)
    s2w <- max(stats::var(ew), 1e-4)
    s2b <- max(stats::var(am), 1e-4)
    d0  <- if (intercept_zu) c(log(s2w / 2), rep(0, pu - 1L)) else rep(0.1, pu)
    g0  <- if (intercept_zh) c(log(s2b / 2), rep(0, qh - 1L)) else rep(0.1, qh)
    start_v <- unname(c(sqrt(s2w / 2), sqrt(s2b / 2), b0, d0, g0))
  }

  R <- if (is.numeric(halton_num)) as.integer(halton_num) else ceiling(sqrt(n)) + 100L
  if (draws == "shared") {
    D  <- .make_draws(R, discard = halton_discard, seed = rand.gtre)
    zr <- D[, 1]; zh <- D[, 2]
  } else {
    D  <- .make_draws_perfirm(N, R, seed = rand.gtre)
    zr <- D$r; zh <- D$h
  }

  i_b <- 3L:(2L + kx)
  i_d <- (3L + kx):(2L + kx + pu)
  i_g <- (3L + kx + pu):(2L + kx + pu + qh)
  CLIP <- .zh_constants$LIN_CLIP

  obj <- .zh_make_objective(Yv, X, Zu, W, gid = match(id, indiv),
                            zr = zr, zh = zh, s = s)
  fn  <- function(p) obj(p)$value

  last <- new.env()
  fn_g <- function(p) {
    r <- obj(p, deriv = TRUE)
    last$p <- p; last$g <- r$grad
    r$value
  }
  gr <- function(p) {
    if (identical(p, last$p)) return(last$g)
    obj(p, deriv = TRUE)$grad
  }

  lower <- c(rep(.zh_constants$SIG_LOWER, 2L), rep(-Inf, npar - 2L))
  p <- start_v
  if (verbose) cat("initial nll:", fn(p), "\n")

  if (maxit.bobyqa > 0 && requireNamespace("minqa", quietly = TRUE)) {
    bo <- try(minqa::bobyqa(p, fn, lower = lower,
                            control = list(maxfun = max(maxit.bobyqa,
                                                        10L * npar^2 + 50L))),
              silent = TRUE)
    if (!inherits(bo, "try-error") && is.finite(bo$fval) && bo$fval < fn(p)) p <- bo$par
    if (verbose) cat("after bobyqa, nll:", fn(p), "\n")
  }

  if (isTRUE(PSopt) && requireNamespace("pso", quietly = TRUE)) {
    if (!is.null(rand.psoptim)) set.seed(rand.psoptim)
    ps <- try(pso::psoptim(p, fn,
                           lower   = pmax(lower, p - 10),
                           upper   = p + 10,
                           control = list(maxit = maxit.psoptim)), silent = TRUE)
    if (!inherits(ps, "try-error") && is.finite(ps$value) && ps$value < fn(p)) p <- ps$par
    if (verbose) cat("after psoptim, nll:", fn(p), "\n")
  }

  f_opt <- if (isTRUE(analytic_grad)) fn_g else fn
  g_opt <- if (isTRUE(analytic_grad)) gr else NULL
  op <- if (Method == "L-BFGS-B") {
    stats::optim(p, f_opt, g_opt, method = Method, lower = lower,
                 hessian = optHessian,
                 control = list(maxit = maxit.optim, trace = if (verbose) 1 else 0))
  } else {
    stats::optim(p, f_opt, g_opt, method = Method, hessian = optHessian,
                 control = list(maxit = maxit.optim, trace = if (verbose) 1 else 0))
  }
  if (verbose) cat("final nll:", op$value, " convergence:", op$convergence, "\n")

  grad_opt <- obj(op$par, deriv = TRUE)$grad
  names(grad_opt) <- pnames
  hess_min_eig <- if (isTRUE(optHessian) && !is.null(op$hessian) &&
                      all(is.finite(op$hessian))) {
    min(eigen((op$hessian + t(op$hessian)) / 2, symmetric = TRUE,
              only.values = TRUE)$values)
  } else NA_real_

  se <- rep(NA_real_, npar)
  if (isTRUE(optHessian) && !is.null(op$hessian) && all(is.finite(op$hessian))) {
    Vc <- try(solve(op$hessian), silent = TRUE)
    if (!inherits(Vc, "try-error")) {
      dv <- diag(Vc)
      se <- ifelse(dv > 0, sqrt(dv), NA_real_)
    }
  }
  out <- rbind(par = op$par, st_err = se, `t-val` = op$par / se)
  colnames(out) <- pnames

  est  <- op$par
  sigv <- est[1]; sigr <- est[2]
  beta <- est[i_b]; delta <- est[i_d]; gam <- est[i_g]
  su_all <- sqrt(exp(pmin(pmax(as.numeric(Zu %*% delta), -CLIP), CLIP)))
  sh_vec <- sqrt(exp(pmin(pmax(as.numeric(W  %*% gam),   -CLIP), CLIP)))
  names(sh_vec) <- as.character(indiv)

  U <- H <- O <- NULL
  if (isTRUE(compute_te)) {
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      warning("mvtnorm not installed: efficiency scores not computed.")
    } else {
      e_all <- s * (Yv - as.numeric(X %*% beta))
      te <- .zh_efficiency(e_all, idx, sigv, sigr, sh_vec, su_all, verbose)
      U <- te$U; H <- te$H; O <- te$O
    }
  }

  results <- list(out          = t(out),
                  opt          = op,
                  data         = data,
                  total_time   = difftime(Sys.time(), t0, units = "secs"),
                  start_v      = start_v,
                  model_name   = "GTRE_ZH",
                  formula      = formula,
                  inefdec_n    = s,
                  U            = U,
                  H            = H,
                  O            = O,
                  sigma_u      = su_all,
                  sigma_h      = sh_vec,
                  loglik       = -op$value,
                  grad         = grad_opt,
                  max_abs_grad = max(abs(grad_opt)),
                  hess_min_eig = hess_min_eig,
                  n = n, N = N, t_i = t_i, R_draws = R,
                  draws        = draws,
                  individual   = individual,
                  index        = idx,
                  zu_names     = colnames(Zu),
                  zh_names     = colnames(W),
                  coefficients = out["par", ],
                  std.errors   = out["st_err", ],
                  t.values     = out["t-val", ],
                  call         = call)
  class(results) <- c("psfm_gtre_zh", "sfareg")
  results
}

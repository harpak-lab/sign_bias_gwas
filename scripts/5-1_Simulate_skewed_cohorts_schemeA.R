############################################################
## Simulation of skewed cohorts, scheme A
##
## A polygenic trait with balanced positive and negative allelic
## effects is simulated in the full population. Cohorts are then
## sampled from increasingly right-skewed portions of the phenotype
## distribution.
##
## Parameters:
##   N_POP        population size
##   M_PAIRS      number of matched +/- SNP pairs
##   H2           target SNP heritability
##   N_SAMP       cohort sample size
##   GAMMA_FIXED  concentration near truncation threshold q
##   TAU_FIXED    mixture weight for uniform bin sampling
##   MAF_MAX      MAF threshold for evaluation
##   BLOCK_SIZE   pseudo-block size for min-p winner selection
##   SKEW_TOL     tolerance for target skew matching
##   MAX_TRIES    maximum sampling attempts per replicate
############################################################

# Libraries
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ashr)
library(cowplot)
library(splines)

## Parameters, input
N_POP        <- 2000000
M_PAIRS      <- 2000
H2           <- 0.6
N_BINS       <- 200

TARGETS      <- 1:10
N_REP        <- 20
N_SAMP       <- 10000

TAU_FIXED    <- 0.01
GAMMA_FIXED  <- 20
Q_GRID       <- seq(0, 0.95, by = 0.0001)

TUNE_SEED    <- 23276
BASE_SEED    <- 1000
GWAS_CHUNK   <- 200
BLOCK_SIZE   <- 6L
CI_LEVEL     <- 0.95

MAF_MAX      <- 0.01
SKEW_TOL     <- 0.5
MAX_TRIES    <- 50

### Calculate skewness (standardized 3rd moment)
skew3 <- function(x){
  m <- mean(x); s <- sd(x)
  if (s == 0) 0 else mean(((x - m) / s)^3)
}

### Build phenotype w/ corresponding genotypes
build_population <- function(N_pop = 140000,
                             M_pairs = 1000,
                             maf_min = 5e-4,
                             maf_max = 8e-3,
                             h2 = 0.6,
                             # demo/focal pair
                             add_demo = TRUE,
                             p_demo = 0.30,
                             beta_demo = 0.40,
                             chunk = 200,
                             seed = 1) {
  set.seed(seed)
  
  M_bg   <- 2 * M_pairs
  M_demo <- if (isTRUE(add_demo)) 2L else 0
  M      <- M_bg + M_demo
  
  # paired MAFs
  p_pair <- 10^(runif(M_pairs, log10(maf_min), log10(maf_max)))
  p_bg   <- rep(p_pair, each = 2)
  
  # demo genetic variance under raw dosages:
  # demo_effect = beta_demo * (G_inc - G_dec)
  # Var(G_inc - G_dec) = 2*Var(G) = 4 p (1-p)
  if (isTRUE(add_demo)) {
    Vg_demo <- (beta_demo^2) * (4 * p_demo * (1 - p_demo))
    if (Vg_demo >= h2) {
      stop(sprintf("beta_demo too large for target h2 (Vg_demo=%.3f >= h2=%.3f). Lower beta_demo or raise h2.",
                   Vg_demo, h2))
    }
  } else {
    Vg_demo <- 0
  }
  h2_bg <- h2 - Vg_demo
  
  # choose ONE per-allele background effect magnitude so that
  # Var(sum_j beta_j*G_j) ≈ h2_bg (independent SNPs; Var(G)=2p(1-p))
  sum_varG_bg <- sum(2 * p_bg * (1 - p_bg))
  bmag_bg <- sqrt(h2_bg / sum_varG_bg)
  
  # balanced +/- effects
  beta_bg <- rep(c(+bmag_bg, -bmag_bg), times = M_pairs)
  
  # full p and beta vectors
  if (isTRUE(add_demo)) {
    p    <- c(p_bg, p_demo, p_demo)
    beta <- c(beta_bg, +beta_demo, -beta_demo)
    demo_inc_idx <- M_bg + 1L
    demo_dec_idx <- M_bg + 2L
  } else {
    p    <- p_bg
    beta <- beta_bg
    demo_inc_idx <- NA
    demo_dec_idx <- NA
  }
  
  # store raw genotypes
  G_raw <- matrix(as.raw(0), nrow = N_pop, ncol = M)
  
  # accumulate genetic value using RAW dosages (NO centering)
  gval <- numeric(N_pop)
  
  for (j0 in seq(1, M_bg, by = chunk)) {
    j1 <- min(M_bg, j0 + chunk - 1)
    jj <- j0:j1
    m  <- length(jj)
    
    probs <- rep(p_bg[jj], each = N_pop)
    Gi <- matrix(rbinom(N_pop * m, 2, probs), nrow = N_pop, ncol = m)
    
    G_raw[, jj] <- matrix(as.raw(Gi), nrow = N_pop, ncol = m)
    
    gval <- gval + as.vector(Gi %*% beta_bg[jj])
  }
  
  # add demo pair
  if (isTRUE(add_demo)) {
    G_demo_inc <- rbinom(N_pop, 2, p_demo)
    G_demo_dec <- rbinom(N_pop, 2, p_demo)
    
    G_raw[, demo_inc_idx] <- as.raw(G_demo_inc)
    G_raw[, demo_dec_idx] <- as.raw(G_demo_dec)
    
    demo_effect <- beta_demo * (G_demo_inc - G_demo_dec)
    gval <- gval + demo_effect
  } else {
    G_demo_inc <- NULL
    G_demo_dec <- NULL
  }
  
  # add noise and build phenotype (NO standardization)
  eps   <- rnorm(N_pop, sd = sqrt(1 - h2))
  Y_pop <- as.numeric(gval + eps)
  
  list(
    N_pop = N_pop,
    M = M,
    M_bg = M_bg,
    p = p,
    beta = beta,
    G_raw = G_raw,
    Y_pop = Y_pop,
    demo_inc_idx = demo_inc_idx,
    demo_dec_idx = demo_dec_idx,
    G_demo_inc = G_demo_inc,
    G_demo_dec = G_demo_dec,
    h2 = h2,
    h2_bg = h2_bg,
    Vg_demo = Vg_demo,
    bmag_bg = bmag_bg
  )
}

# Precompute quantile bins for pop, fast cohort sampling
prep_bins <- function(Y, B = 200) {
  probs <- seq(0, 1, length.out = B + 1)
  brks <- as.numeric(quantile(Y, probs))
  brks <- unique(brks)
  B_eff <- length(brks) - 1
  if (B_eff < 10) stop("Too few unique quantile breaks; reduce B.")
  
  bin <- cut(Y, breaks = brks, include.lowest = TRUE, labels = FALSE)
  idx_by_bin <- split(seq_along(Y), bin)
  sizes <- vapply(idx_by_bin, length, integer(1))
  r_mid <- (seq_len(B_eff) - 0.5) / B_eff
  
  list(Y = Y, brks = brks, idx_by_bin = idx_by_bin, sizes = sizes, r_mid = r_mid, B = B_eff)
}

allocate_counts <- function(N_samp, prob, sizes) {
  prob <- prob / sum(prob)
  cnt <- as.vector(rmultinom(1, N_samp, prob = prob))
  
  overflow <- which(cnt > sizes)
  while (length(overflow) > 0) {
    excess <- sum(cnt[overflow] - sizes[overflow])
    cnt[overflow] <- sizes[overflow]
    
    cap <- sizes - cnt
    if (sum(cap) == 0) stop("No remaining capacity to allocate sample; reduce N_samp or B.")
    p2 <- cap * prob
    p2 <- p2 / sum(p2)
    add <- as.vector(rmultinom(1, excess, prob = p2))
    cnt <- cnt + add
    
    overflow <- which(cnt > sizes)
  }
  cnt
}

### Cohort sampler 
sample_cohort <- function(prep,
                          N_samp = 7000,
                          q = 0,
                          gamma = 5,
                          tau = 0.05,
                          seed = 1) {
  set.seed(seed)
  
  r <- prep$r_mid
  keep_bins <- which(r >= q)
  
  if (length(keep_bins) < 5) stop("q too high")
  avail <- sum(prep$sizes[keep_bins])
  if (N_samp > avail) stop(sprintf("not enough inds above q: avail=%d, need=%d", avail, N_samp))
  
  r_keep <- (r[keep_bins] - q) / (1 - q + 1e-12)
  
  # weights high near cutoff, low near top
  w <- (1 - r_keep)^gamma
  p_keep <- w / sum(w)
  
  K <- length(keep_bins)
  prob <- rep(0, prep$B)
  prob[keep_bins] <-tau * rep(1 / K, K) + (1 -tau) * p_keep
  
  cnt <- allocate_counts(N_samp, prob, prep$sizes)
  
  idx <- unlist(mapply(function(v, k) if (k > 0) sample(v, k, FALSE) else integer(0),
                       prep$idx_by_bin, cnt,
                       SIMPLIFY = FALSE, USE.NAMES = FALSE))
  y_samp <- prep$Y[idx]
  
  list(idx = idx, Y_samp = y_samp,
       mean = mean(y_samp), skew = skew3(y_samp),
       q = q, gamma = gamma,tau =tau)
}

### GWAS (linreg)
gwas_all_snps <- function(pop, idx, chunk = 200, ...) {
  stopifnot(!is.null(pop$G_raw), !is.null(pop$Y_pop))
  n <- length(idx)
  y <- pop$Y_pop[idx]
  y <- y - mean(y)
  Syy <- sum(y^2)
  
  M <- ncol(pop$G_raw)
  
  betahat <- rep(NA, M)
  se      <- rep(NA, M)
  maf     <- rep(NA, M)
  mac     <- rep(NA, M)
  af      <- rep(NA, M)
  
  for (j0 in seq(1, M, by = chunk)) {
    j1 <- min(M, j0 + chunk - 1)
    jj <- j0:j1
    m  <- length(jj)
    
    G <- matrix(as.integer(pop$G_raw[idx, jj, drop = FALSE]), nrow = n, ncol = m)
    
    gsum <- colSums(G)
    af_j <- gsum / (2 * n)
    maf_j <- pmin(af_j, 1 - af_j)
    mac_j <- pmin(gsum, 2*n - gsum)
    
    af[jj]  <- af_j
    maf[jj] <- maf_j
    mac[jj] <- mac_j
    
    keep <- (mac_j >= 1) & (maf_j > 0)
    if (!any(keep)) next
    
    gbar <- gsum / n
    Sxx  <- colSums(G^2) - n * (gbar^2)
    Sxy  <- as.numeric(crossprod(G, y))
    
    b <- Sxy / Sxx
    
    # orient effect to minor allele
    flip <- af_j > 0.5
    b[flip] <- -b[flip]
    
    SSE <- Syy - (Sxy^2) / Sxx
    sigma2 <- SSE / (n - 2)
    s <- sqrt(sigma2 / Sxx)
    
    betahat[jj] <- b
    se[jj]      <- s
  }
  
  z <- betahat / se
  data.frame(snp = seq_len(M),
             af = af,
             maf = maf,
             mac = mac,
             betahat = betahat,
             se = se,
             z = z)
}

# ASH per-SNP sign bias: sb = Pr(beta>0) - Pr(beta<0)
# =========================================================
run_ash_sb <- function(bhat, se) {
  if (!requireNamespace("ashr", quietly = TRUE)) {
    stop("Please install ashr: install.packages('ashr')")
  }
  sb <- rep(NA, length(bhat))
  ok <- is.finite(bhat) & is.finite(se) & (se > 0)
  if (!any(ok)) return(sb)
  
  sb_ok <- tryCatch({
    fit <- ashr::ash(betahat = bhat[ok], sebetahat = se[ok],
                     method = "fdr",
                     mixcompdist = "normal",
                     optmethod = "mixSQP")
    
    if (!is.null(fit$result) &&
        all(c("PositiveProb", "NegativeProb") %in% colnames(fit$result))) {
      as.numeric(fit$result$PositiveProb - fit$result$NegativeProb)
    } else {
      p_ge0 <- as.numeric(ashr::get_posterior_prob(fit, l = 0,     u = Inf))
      p_le0 <- as.numeric(ashr::get_posterior_prob(fit, l = -Inf,  u = 0))
      as.numeric(p_ge0 - p_le0)  # == Pr(beta>0) - Pr(beta<0)
    }
  }, error = function(e) {
    z <- bhat[ok] / se[ok]
    sign(z) * (2 * pnorm(abs(z)) - 1)
  })
  
  sb[ok] <- sb_ok
  sb
}

# Cohort-level sign bias from per-SNP sb values:
# S = sum(sb)/sum(|sb|)
cohort_sign_bias_from_sb <- function(sb_vec) {
  ok <- is.finite(sb_vec)
  if (!any(ok)) return(NA)
  den <- sum(abs(sb_vec[ok]), na.rm = TRUE)
  if (!is.finite(den) || den <= .Machine$double.eps) return(NA)
  sum(sb_vec[ok], na.rm = TRUE) / den
}

# TRUE sign for a SNP, minor-allele oriented (unchanged)
true_minor_sign <- function(af, beta_true) {
  s <- sign(beta_true)
  ifelse(af <= 0.5, s, -s)
}

.safe_ratio <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || den <= .Machine$double.eps) return(NA)
  num / den
}

# Diagnostics, computed on SAME MAF subset (<= maf_max) and MAC>=mac_min
gwas_diag_ratios <- function(gwas, mac_min = 1, maf_max = Inf) {
  ok <- is.finite(gwas$betahat) & is.finite(gwas$se) & (gwas$se > 0) &
    is.finite(gwas$z) & is.finite(gwas$mac) & (gwas$mac >= mac_min) &
    is.finite(gwas$maf) & (gwas$maf > 0) & (gwas$maf <= maf_max)
  
  if (!any(ok)) {
    return(list(
      ratio_absb = NA, ratio_se = NA, ratio_t2 = NA,
      n_inc = 0, n_dec = 0
    ))
  }
  
  gg <- gwas[ok, , drop = FALSE]
  
  inc <- gg$betahat > 0
  dec <- gg$betahat < 0
  
  n_inc <- sum(inc, na.rm = TRUE)
  n_dec <- sum(dec, na.rm = TRUE)
  
  if (n_inc == 0 || n_dec == 0) {
    return(list(
      ratio_absb = NA, ratio_se = NA, ratio_t2 = NA,
      n_inc = n_inc, n_dec = n_dec
    ))
  }
  
  absb_inc <- mean(abs(gg$betahat[inc]), na.rm = TRUE)
  absb_dec <- mean(abs(gg$betahat[dec]), na.rm = TRUE)
  
  se_inc <- mean(gg$se[inc], na.rm = TRUE)
  se_dec <- mean(gg$se[dec], na.rm = TRUE)
  
  t2 <- gg$z^2
  t2_inc <- mean(t2[inc], na.rm = TRUE)
  t2_dec <- mean(t2[dec], na.rm = TRUE)
  
  list(
    ratio_absb = .safe_ratio(absb_inc, absb_dec),
    ratio_se   = .safe_ratio(se_inc,   se_dec),
    ratio_t2   = .safe_ratio(t2_inc,   t2_dec),
    n_inc = n_inc,
    n_dec = n_dec
  )
}

# =========================================================
# pseudo-block min-p selection among a pool of SNP indices
#   - blocks defined by SNP index: block_id = floor((snp-1)/block_size)+1
#   - select the SNP with smallest GWAS p (two-sided) within each block
# =========================================================
.minp_per_block <- function(gwas, idx_pool, block_size = 6L) {
  if (length(idx_pool) == 0) return(integer(0))
  
  p <- rep(Inf, length(idx_pool))
  zsub <- gwas$z[idx_pool]
  okp <- is.finite(zsub)
  p[okp] <- 2 * stats::pnorm(-abs(zsub[okp]))
  
  block_id <- ((gwas$snp[idx_pool] - 1L) %/% as.integer(block_size)) + 1L
  spl <- split(seq_along(idx_pool), block_id)
  
  winners_local <- vapply(spl, function(ii_local) {
    ii_local[which.min(p[ii_local])]
  }, integer(1))
  
  idx_pool[winners_local]
}


# sample until skew within tolerance 
# =========================================================
.pick_sample_with_tol <- function(prep, N_samp, q, gamma,tau, seed,
                                  target_sk,
                                  skew_tol = Inf,
                                  max_tries = 1,
                                  give_sample = NULL) {
  best <- NULL
  best_err <- Inf
  
  if (!is.null(give_sample)) {
    sk0 <- skew3(give_sample$Y_samp)
    err0 <- abs(sk0 - target_sk)
    best <- give_sample
    best_err <- err0
    if (is.finite(err0) && err0 <= skew_tol) {
      return(list(samp = give_sample, skew = sk0, err = err0, in_tol = TRUE, n_tries = 1L))
    }
  }
  
  start_try <- if (!is.null(give_sample)) 2L else 1L
  if (max_tries < start_try) {
    if (is.null(best)) return(list(samp = NULL, skew = NA, err = NA, in_tol = FALSE, n_tries = max_tries))
    return(list(samp = best, skew = skew3(best$Y_samp), err = best_err, in_tol = FALSE, n_tries = max_tries))
  }
  
  for (tt in start_try:max_tries) {
    s <- tryCatch(
      sample_cohort(prep,
                    N_samp = N_samp,
                    q = q,
                    gamma = gamma,
                    tau =tau,
                    seed = seed + 1000000 * (tt - 1L)),
      error = function(e) NULL
    )
    if (is.null(s)) next
    sk <- skew3(s$Y_samp)
    err <- abs(sk - target_sk)
    
    if (is.finite(err) && err < best_err) {
      best <- s
      best_err <- err
    }
    
    if (is.finite(err) && err <= skew_tol) {
      return(list(samp = s, skew = sk, err = err, in_tol = TRUE, n_tries = tt))
    }
  }
  
  if (is.null(best)) return(list(samp = NULL, skew = NA, err = NA, in_tol = FALSE, n_tries = max_tries))
  list(samp = best, skew = skew3(best$Y_samp), err = best_err, in_tol = FALSE, n_tries = max_tries)
}

# Fit ASH on ALL SNPs with MAC>=1
one_cohort_bias <- function(pop, prep,
                            target_sk, q, gamma,
                            N_samp,
                            tau = 0.05,
                            seed = 1,
                            gwas_chunk = 200,
                            mac_min = 0,
                            select = "block_min",
                            block_size = 6,
                            n_top = 1000,
                            maf_max = Inf,
                            give_sample = NULL,
                            skew_tol = Inf,
                            max_tries = 1) {
  
  pick <- .pick_sample_with_tol(
    prep = prep,
    N_samp = N_samp,
    q = q,
    gamma = gamma,
    tau =tau,
    seed = seed,
    target_sk = target_sk,
    skew_tol = skew_tol,
    max_tries = max_tries,
    give_sample = give_sample
  )
  
  if (is.null(pick$samp)) {
    return(tibble::tibble(
      target_sk = target_sk,
      skew_obs  = NA,
      skew_err  = NA,
      in_tol    = FALSE,
      n_tries   = max_tries,
      q         = q,
      gamma     = gamma,
      bias_sig  = NA,
      true_bias_sample = NA,
      n_used    = 0,
      seed      = seed,
      ratio_absb = NA,
      ratio_se   = NA,
      ratio_t2   = NA,
      n_inc_gwas = 0,
      n_dec_gwas = 0
    ))
  }
  
  samp <- pick$samp
  
  if (!isTRUE(pick$in_tol)) {
    return(tibble::tibble(
      target_sk = target_sk,
      skew_obs  = pick$skew,
      skew_err  = pick$err,
      in_tol    = FALSE,
      n_tries   = pick$n_tries,
      q         = q,
      gamma     = gamma,
      bias_sig  = NA,
      true_bias_sample = NA,
      n_used    = 0,
      seed      = seed,
      ratio_absb = NA,
      ratio_se   = NA,
      ratio_t2   = NA,
      n_inc_gwas = 0,
      n_dec_gwas = 0
    ))
  }
  
  gwas <- gwas_all_snps(pop, idx = samp$idx, chunk = gwas_chunk, mac_min = mac_min)
  
  # diagnostics computed on the SAME maf_max subset
  diag <- gwas_diag_ratios(gwas, mac_min = 1, maf_max = maf_max)
  
  # FIT ASH on ALL SNPs with MAC>=1
  ok_fit <- is.finite(gwas$betahat) & is.finite(gwas$se) & (gwas$se > 0) &
    is.finite(gwas$mac) & (gwas$mac >= 1)
  
  sb_full <- rep(NA, nrow(gwas))
  if (any(ok_fit)) {
    sb_full[ok_fit] <- run_ash_sb(gwas$betahat[ok_fit], gwas$se[ok_fit])
  }
  
  # Eligible set: MAF<=maf_max, MAC>=1, sb finite 
  ok_eval <- ok_fit &
    is.finite(gwas$maf) & (gwas$maf > 0) & (gwas$maf <= maf_max) &
    is.finite(sb_full) &
    is.finite(gwas$z)
  
  # TRUE bias
  true_bias_sample <- if (!any(ok_eval)) {
    NA
  } else {
    beta_true <- pop$beta[gwas$snp[ok_eval]]
    s_true_minor <- true_minor_sign(gwas$af[ok_eval], beta_true)
    mean(s_true_minor)
  }
  
  # INFERRED bias
  idx_pool <- which(ok_eval)
  idx_win  <- .minp_per_block(gwas, idx_pool, block_size = 6L)
  
  bias_sig <- if (length(idx_win) == 0) {
    NA
  } else {
    cohort_sign_bias_from_sb(sb_full[idx_win])
  }
  
  tibble::tibble(
    target_sk = target_sk,
    skew_obs  = pick$skew,
    skew_err  = pick$err,
    in_tol    = TRUE,
    n_tries   = pick$n_tries,
    q         = q,
    gamma     = gamma,
    bias_sig  = bias_sig,
    true_bias_sample = true_bias_sample,
    n_used    = length(idx_win),   # number of winners used for inferred bias
    seed      = seed,
    ratio_absb = diag$ratio_absb,
    ratio_se   = diag$ratio_se,
    ratio_t2   = diag$ratio_t2,
    n_inc_gwas = diag$n_inc,
    n_dec_gwas = diag$n_dec
  )
}

## Tune ONLY q (gamma +tau fixed)

tune_q_for_skew <- function(prep,
                            target_sk,
                            N_samp,
                            gamma_fixed,
                            tau_fixed,
                            q_grid = seq(0, 0.95, by = 0.01),
                            seed = 1,
                            skew_tol = Inf) {
  best <- NULL
  best_err <- Inf
  
  for (q in q_grid) {
    samp <- tryCatch(
      sample_cohort(prep,
                    N_samp = N_samp,
                    q = q,
                    gamma = gamma_fixed,
                    tau =tau_fixed,
                    seed = seed),
      error = function(e) NULL
    )
    if (is.null(samp)) next
    
    sk <- skew3(samp$Y_samp)
    err <- abs(sk - target_sk)
    
    if (is.finite(err) && err < best_err) {
      best_err <- err
      best <- list(q = q, skew = sk)
    }
  }
  
  if (is.null(best)) {
    return(tibble(target_sk = target_sk,
                  q = NA,
                  gamma = gamma_fixed,
                  tau =tau_fixed,
                  skew_pilot = NA,
                  err = NA,
                  in_tol = FALSE))
  }
  
  in_tol <- is.finite(best_err) && (best_err <= skew_tol)
  
  tibble(
    target_sk = target_sk,
    q = if (in_tol) best$q else NA,
    gamma = gamma_fixed,
    tau =tau_fixed,
    skew_pilot = best$skew,
    err = best_err,
    in_tol = in_tol
  )
}


## simulate with q-only tuning (one point per target)
simulate_signbias_vs_skew_qonly <- function(pop, prep,
                                            targets = 1:10,
                                            n_rep = 20,
                                            N_samp = 10000,
                                            tau_fixed = 0.01,
                                            gamma_fixed = 20,
                                            tune_seed = 1,
                                            q_grid = seq(0, 0.95, by = 0.01),
                                            gwas_chunk = 200,
                                            mac_min = 0,
                                            select = "block_min",
                                            block_size = 6,
                                            n_top = 1000,
                                            maf_max = Inf,
                                            ci_level = 0.95,
                                            base_seed = 1000,
                                            skew_tol = Inf,
                                            max_tries = 1) {
  
  tuning <- bind_rows(lapply(seq_along(targets), function(i) {
    tune_q_for_skew(prep,
                    target_sk = targets[i],
                    N_samp = N_samp,
                    gamma_fixed = gamma_fixed,
                    tau_fixed =tau_fixed,
                    q_grid = q_grid,
                    seed = tune_seed + i,
                    skew_tol = skew_tol)
  }))
  
  reps <- bind_rows(lapply(seq_along(targets), function(i) {
    tsk <- targets[i]
    row <- tuning %>% filter(target_sk == tsk)
    
    if (!is.finite(row$q[1])) {
      return(tibble(
        target_sk = tsk,
        skew_obs = NA,
        skew_err = NA,
        in_tol = FALSE,
        n_tries = max_tries,
        q = NA,
        gamma = gamma_fixed,
        bias_sig = NA,
        true_bias_sample = NA,
        n_used = 0,
        seed = NA,
        ratio_absb = NA,
        ratio_se = NA,
        ratio_t2 = NA,
        n_inc_gwas = 0,
        n_dec_gwas = 0
      ))
    }
    
    bind_rows(lapply(seq_len(n_rep), function(r) {
      one_cohort_bias(pop, prep,
                      target_sk = tsk,
                      q = row$q[1],
                      gamma = gamma_fixed,
                      N_samp = N_samp,
                      tau =tau_fixed,
                      seed = base_seed + 10000 * i + r,
                      gwas_chunk = gwas_chunk,
                      mac_min = mac_min,
                      select = select,
                      block_size = block_size,
                      n_top = n_top,
                      maf_max = maf_max,
                      give_sample = NULL,
                      skew_tol = skew_tol,
                      max_tries = max_tries)
    }))
  }))
  
  alpha <- 1 - ci_level
  summ <- reps %>%
    group_by(target_sk) %>%
    summarise(
      n_rep = sum(is.finite(bias_sig)),
      n_in_tol = sum(isTRUE(in_tol)),
      mean_skew_obs = mean(skew_obs, na.rm = TRUE),
      
      mean_bias = mean(bias_sig, na.rm = TRUE),
      sd_bias   = sd(bias_sig, na.rm = TRUE),
      se_bias   = sd_bias / sqrt(pmax(n_rep, 1)),
      
      mean_true = mean(true_bias_sample, na.rm = TRUE),
      sd_true   = sd(true_bias_sample, na.rm = TRUE),
      se_true   = sd_true / sqrt(pmax(n_rep, 1)),
      
      tcrit     = ifelse(n_rep > 1, qt(1 - alpha/2, df = n_rep - 1), NA),
      
      ci_low_bias  = mean_bias - tcrit * se_bias,
      ci_high_bias = mean_bias + tcrit * se_bias,
      
      ci_low_true  = mean_true - tcrit * se_true,
      ci_high_true = mean_true + tcrit * se_true,
      
      .groups = "drop"
    ) %>%
    left_join(tuning, by = "target_sk")
  
  pE_plot <- bind_rows(
    summ %>% transmute(x = target_sk,
                       series = "Cohort true sign bias",
                       mean = mean_true,
                       ymin = ci_low_true,
                       ymax = ci_high_true,
                       shape = "true"),
    summ %>% transmute(x = target_sk,
                       series = "Estimated sign bias",
                       mean = mean_bias,
                       ymin = ci_low_bias,
                       ymax = ci_high_bias,
                       shape = "est")
  )
  
  pE <- ggplot(pE_plot, aes(x = x, y = mean, color = series, group = series)) +
    geom_hline(yintercept = 0, linewidth = 1, color = "#D5B60A") +
    geom_line(linewidth = 0.7) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.035, linewidth = 0.5) +
    geom_point(aes(shape = shape), size = 3.0, stroke = 0.9) +
    scale_shape_manual(values = c(true = 17, est = 16), guide = "none") +
    scale_color_manual(values = c(
      "Cohort true sign bias" = "purple4",
      "Estimated sign bias"   = "black"
    )) +
    scale_x_continuous(breaks = pretty(summ$target_sk), limits = range(summ$target_sk)) +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_text(size = 12),
      legend.position = "none"
    )
  
  list(tuning = tuning, reps = reps, summary = summ, plot = pE)
}

############################################################
## RUN
############################################################
pop <- build_population(
  N_pop = N_POP,
  M_pairs = M_PAIRS,
  maf_max = 0.1,
  add_demo = FALSE,
  seed = 1
)

prep <- prep_bins(pop$Y_pop, B = N_BINS)

targets <- TARGETS

out <- simulate_signbias_vs_skew_qonly(
  pop, prep,
  targets = targets,
  n_rep = N_REP,
  N_samp = N_SAMP,
  tau_fixed = TAU_FIXED,
  gamma_fixed = GAMMA_FIXED,
  q_grid = Q_GRID,
  tune_seed = TUNE_SEED,
  maf_max = MAF_MAX,
  skew_tol = SKEW_TOL,
  max_tries = MAX_TRIES
)

print(out$tuning)
print(out$summary)
print(out$plot)

############################################################
#baseline skew=0 reps + diagnostic pDiag (for specified MAF_MAX)
############################################################

reps0 <- dplyr::bind_rows(lapply(seq_len(N_REP), function(r) {
  samp0 <- sample_cohort(prep, N_samp = N_SAMP, q = 0, gamma = 0,tau = 0.01, seed = 900000 + r)
  one_cohort_bias(pop, prep,
                  target_sk = 0,
                  q = 0, gamma = 0,
                  N_samp = N_SAMP,
                  tau = 1,
                  seed = 900000 + r,
                  maf_max = MAF_MAX,
                  give_sample = samp0,
                  skew_tol = SKEW_TOL,
                  max_tries = MAX_TRIES)
}))

reps_all <- dplyr::bind_rows(reps0, out$reps)

ci_level <- CI_LEVEL
alpha <- 1 - ci_level

diag_summ <- reps_all %>%
  group_by(target_sk) %>%
  summarise(
    n_rep = n(),
    
    mean_absb = mean(ratio_absb, na.rm = TRUE),
    sd_absb   = sd(ratio_absb, na.rm = TRUE),
    se_absb   = sd_absb / sqrt(sum(is.finite(ratio_absb))),
    
    mean_se = mean(ratio_se, na.rm = TRUE),
    sd_se   = sd(ratio_se, na.rm = TRUE),
    se_se   = sd_se / sqrt(sum(is.finite(ratio_se))),
    
    mean_t2 = mean(ratio_t2, na.rm = TRUE),
    sd_t2   = sd(ratio_t2, na.rm = TRUE),
    se_t2   = sd_t2 / sqrt(sum(is.finite(ratio_t2))),
    
    tcrit = ifelse(n_rep > 1, qt(1 - alpha/2, df = n_rep - 1), NA),
    
    lo_absb = mean_absb - tcrit * se_absb,
    hi_absb = mean_absb + tcrit * se_absb,
    
    lo_se = mean_se - tcrit * se_se,
    hi_se = mean_se + tcrit * se_se,
    
    lo_t2 = mean_t2 - tcrit * se_t2,
    hi_t2 = mean_t2 + tcrit * se_t2,
    
    .groups = "drop"
  )

diag_plot_df <- bind_rows(
  diag_summ %>% transmute(target_sk, metric = "|beta| mean ratio (inc/dec)",
                          mean = mean_absb, lo = lo_absb, hi = hi_absb),
  diag_summ %>% transmute(target_sk, metric = "SE mean ratio (inc/dec)",
                          mean = mean_se, lo = lo_se, hi = hi_se),
  diag_summ %>% transmute(target_sk, metric = "T^2 mean ratio (inc/dec)",
                          mean = mean_t2, lo = lo_t2, hi = hi_t2)
) %>%
  mutate(
    target_sk_f = factor(target_sk, levels = rev(sort(unique(target_sk))))
  )

pDiag <- ggplot(diag_plot_df, aes(x = mean, y = target_sk_f)) +
  geom_vline(xintercept = 1, linewidth = 0.8, color = "#D5B60A") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.6) +
  geom_point(size = 2.2) +
  facet_wrap(~ metric, ncol = 1, scales = "free_x") +
  theme_classic(base_size = 12) +
  labs(
    x = "Ratio (increasing / decreasing)",
    y = "Target skew",
    title = paste0("GWAS diagnostic ratios computed within MAF ≤ ", MAF_MAX)
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 11)
  )


#####

demo_cohort <- sample_cohort(
  prep,
  N_samp = N_SAMP,
  q      = 0.46,
  gamma  = GAMMA_FIXED,
  tau    = TAU_FIXED,
  seed   = 199
)
demo_cohort$skew


##### Plotting

## Pop
dens_pop2 <- density(pop$Y_pop)
dd_pop2 <- data.frame(x=dens_pop2$x, 
                      y=dens_pop2$y)

pA <- ggplot(dd_pop2, aes(x,y)) + 
  geom_area(fill="#FDE992", color=NA) +
  geom_line(color="black", linewidth=0.6) +
  geom_vline(
    aes(xintercept = -0.1),
    linewidth = 0.8,
    linetype = 2
  ) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10))

## Cohort
Y_demo2 <- demo_cohort$Y_samp
dens_samp2 <- density(Y_demo2, adjust=3)
dd_samp2 <- data.frame(x=dens_samp2$x, 
                       y=dens_samp2$y)

pB<- ggplot(dd_samp2, aes(x,y)) + 
  geom_area(fill="#EFEFEF", color=NA) +
  geom_line(color="black", linewidth=0.6) + 
  theme_classic() +
  theme(axis.title.y     = element_blank(),
        axis.title.x     = element_blank(),
        axis.text.x      = element_text(size = 10),
        axis.text.y      = element_text(size = 10)) 

## Skew
res_sim <- out$summary

# skew=0
skew0 <- reps_all %>% 
  filter(target_sk == 0) %>%
  summarise(
    target_sk = 0,
    n_rep = sum(is.finite(bias_sig)),
    n_in_tol = sum(isTRUE(in_tol)),
    mean_skew_obs = mean(skew_obs, na.rm = TRUE),
    
    mean_bias = mean(bias_sig, na.rm = TRUE),
    sd_bias   = sd(bias_sig, na.rm = TRUE),
    se_bias   = sd_bias / sqrt(pmax(n_rep, 1)),
    
    mean_true = mean(true_bias_sample, na.rm = TRUE),
    sd_true   = sd(true_bias_sample, na.rm = TRUE),
    se_true   = sd_true / sqrt(pmax(n_rep, 1)),
    
    tcrit     = ifelse(n_rep > 1, qt(1 - alpha/2, df = n_rep - 1), NA),
    
    ci_low_bias  = mean_bias - tcrit * se_bias,
    ci_high_bias = mean_bias + tcrit * se_bias,
    
    ci_low_true  = mean_true - tcrit * se_true,
    ci_high_true = mean_true + tcrit * se_true,
    
    .groups = "drop"
  )

# bind + force column set/order to match res_sim
res_sim <- bind_rows(res_sim, skew0) %>%
  select(all_of(names(res_sim))) %>%
  arrange(target_sk)



pC_plot <- bind_rows(
  res_sim %>% transmute(x = mean_skew_obs,
                        series = "Cohort true sign bias",
                        mean = mean_true,
                        ymin = ci_low_true,
                        ymax = ci_high_true,
                        shape = "true"),
  res_sim %>% transmute(x = mean_skew_obs,
                        series = "Estimated sign bias",
                        mean = mean_bias,
                        ymin = ci_low_bias,
                        ymax = ci_high_bias,
                        shape = "est")
)


pop_pts <- pC_plot %>%
  dplyr::distinct(x) %>%
  dplyr::mutate(
    mean   = 0,
    ymin   = NA,
    ymax   = NA,
    series = "Population sign bias"
  )

pC_plot2 <- dplyr::bind_rows(pC_plot, pop_pts)

p_y <- dplyr::filter(pC_plot2, series == "Population sign bias")
p_p <- dplyr::filter(pC_plot2, series == "Cohort true sign bias")
p_b <- dplyr::filter(pC_plot2, series == "Estimated sign bias")

pC <- ggplot() +
  geom_line(data = p_y, aes(x = x, y = mean, group = 1),
            linewidth = 0.7, color = "#D5B60A") +
  geom_point(data = p_y, aes(x = x, y = mean),
             shape = 16, size = 2, color = "#D5B60A") +
  

  geom_line(data = p_p, aes(x = x, y = mean, group = 1),
            linewidth = 0.7, color = "purple4") +
  geom_errorbar(
    data = dplyr::filter(p_p, !is.na(ymin) & !is.na(ymax)),
    aes(x = x, ymin = ymin, ymax = ymax),
    width = 0.2, linewidth = 0.5, alpha = 0.7, color = "purple4"
  ) +
  geom_point(data = p_p, aes(x = x, y = mean),
             shape = 16, size = 2, color = "purple4") +
  
  geom_line(data = p_b, aes(x = x, y = mean, group = 1),
            linewidth = 0.7, color = "black") +
  geom_errorbar(
    data = dplyr::filter(p_b, !is.na(ymin) & !is.na(ymax)),
    aes(x = x, ymin = ymin, ymax = ymax),
    width = 0.2, linewidth = 0.5, alpha = 0.7, color = "black"
  ) +
  geom_point(data = p_b, aes(x = x, y = mean),
             shape = 16, size = 2, color = "black") +
  
  scale_x_continuous(
    breaks = pretty(out$summary$mean_skew_obs),
    limits = range(out$summary$mean_skew_obs)
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_text(size = 10),
    legend.position = "none") + 
  xlim(1, 10.25) +
  ylim(-0.1, 1)

### plot
top_fig3<-plot_grid(pA,pB,pC, nrow = 1)
top_fig3
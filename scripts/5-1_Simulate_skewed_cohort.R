#!/usr/bin/env Rscript

#####################################################################################
## Simulating skewed phenotypes
#####################################################################################
##
## --single:
##   Runs on a single cohort and outputs results
##     sims_run_f3 <- list(Y_pop, Y_demo, G_inc, 
##   G_dec, gwas_cohort)
##
## --sim:
##   Runs across a specified number of replicates as a 
##   plain text table to --out_sim
##  
#####################################################################################
#
# Modes:
#   --single
#   --sim
# 
# Common:
#   --N_pop (int)        [default 2000000]
#   --pop_seed (int)     [default 1]
#   --chunk (int)        [default 200]
#   --B_bins (int)       [default 200]
#   --out_sample (file)  [default sample.rds]
#   --out_sim (file)     [default out_reps.txt]
# 
# Single cohort:
#   --N_samp (int)       [default 10000]
#   --q (num)            [default 0]
#   --gamma (num)        [default 19.2]
#   --eta (num)          [default 0.1]
#   --sample_seed (int)  [default 199]
#   --panels (bool)      [default false]
#     If true: gwas_cohort is computed from fixed panel GWAS cohort:
#       N_samp=10000, q=0.5, gamma=3.39, eta=0.1, seed=1
# 
# Sim:
#   --n_rep (int)             [default 20]
#   --targets (csv list)      e.g. \"0,0.5,1,1.5\"
#    OR:
#     --targets_start (num)   [default 0]
#     --targets_end (num)     [default 4]
#     --targets_step (num)    [default 0.1]
# 
# Example:
#
#   Panels A-D
#   Rscript simulate_signbias_cli.R --single --N_pop 2000000 --N_samp 10000 --q 0 
#   --gamma 19.2 --eta 0.1 --panels true --sample_seed 199 
#   --out_sample sim1_sample.rds
#
#   Panel E
#   Rscript simulate_signbias_cli.R --sim --N_pop 2000000 --n_rep 20 
#  --targets_start 0 --targets_end 4 --targets_step 0.1 --out_sim sim_1_out_reps.txt
#
#####################################################################################

## CLI arg parsing
.has_flag <- function(args, flag) any(args == flag)

.get_value <- function(args, key, default = NULL) {
  w <- which(args == key)
  if (length(w) == 0) return(default)
  if (w[1] == length(args)) stop("Missing value after ", key)
  args[w[1] + 1]
}

.as_num <- function(x, key) {
  suppressWarnings(v <- as.numeric(x))
  if (!is.finite(v)) stop("Bad numeric for ", key, ": ", x)
  v
}

.as_int <- function(x, key) {
  suppressWarnings(v <- as.integer(x))
  if (is.na(v)) stop("Bad integer for ", key, ": ", x)
  v
}

.as_bool <- function(x, key) {
  x <- tolower(as.character(x))
  if (x %in% c("true","t","1","yes","y")) return(TRUE)
  if (x %in% c("false","f","0","no","n")) return(FALSE)
  stop("Bad boolean for ", key, ": ", x, " (use true/false)")
}

.parse_num_list <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  parts <- strsplit(x, ",", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  v <- suppressWarnings(as.numeric(parts))
  if (any(!is.finite(v))) stop("Could not parse --targets list: ", x)
  v
}

args <- commandArgs(trailingOnly = TRUE)

run_single <- .has_flag(args, "--single")
run_sim    <- .has_flag(args, "--sim")
if (!run_single && !run_sim) run_single <- TRUE

### Calculate skewness (standardized 3rd moment)
skew3 <- function(x){
  m <- mean(x); s <- sd(x)
  if (s == 0) 0 else mean(((x - m) / s)^3)
}


### Build phenotype w/ corresponding genotypes
# N_pop = # individuals in pop
# M_pairs = # of pairs of background SNPs (one + and - in sign)
# maf_min = min MAF
# max_min = max MAF
# h2 = heritability of background SNPs
# add demo = add a pair of focal SNPs for visualization
# p_demo - allele frequency for demonstration SNPs
# beta_demo = effect size magnitude for focal demo SNPs
# p_d = minor allele frequency for demonstration SNPs
# chunk = chunk size (# of SNPs)
# seed = RNG seed value
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
  
  M_bg <- 2 * M_pairs
  M_demo <- if (isTRUE(add_demo)) 2L else 0L
  M <- M_bg + M_demo
  
  p_pair <- 10^(runif(M_pairs, log10(maf_min), log10(maf_max)))
  p_bg   <- rep(p_pair, each = 2)
  
  if (isTRUE(add_demo)) {
    h2_demo <- 2 * beta_demo^2
    if (h2_demo >= h2) {
      stop(sprintf("need lower beta_demo or raise h2.",
                   h2_demo, h2))
    }
  } else {
    h2_demo <- 0
  }
  h2_bg <- h2 - h2_demo
  
  
  # background 
  bmag_bg <- sqrt(h2_bg / M_bg)
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
    demo_inc_idx <- NA_integer_
    demo_dec_idx <- NA_integer_
  }
  
  # store raw genotypes
  G_raw <- matrix(as.raw(0), nrow = N_pop, ncol = M)
  
  # accumulate genetic value: background 
  gval <- numeric(N_pop)
  
  for (j0 in seq(1, M_bg, by = chunk)) {
    j1 <- min(M_bg, j0 + chunk - 1)
    jj <- j0:j1
    m  <- length(jj)
    
    probs <- rep(p_bg[jj], each = N_pop)
    Gi <- matrix(rbinom(N_pop * m, 2, probs), nrow = N_pop, ncol = m)
    
    G_raw[, jj] <- matrix(as.raw(Gi), nrow = N_pop, ncol = m)
    
    mu <- 2 * p_bg[jj]
    sd <- sqrt(2 * p_bg[jj] * (1 - p_bg[jj]))
    Gis <- sweep(Gi, 2, mu, `-`)
    Gis <- sweep(Gis, 2, sd, `/`)
    
    gval <- gval + as.vector(Gis %*% beta_bg[jj])
  }
  
  # add demo pair
  if (isTRUE(add_demo)) {
    mu_d <- 2 * p_demo
    sd_d <- sqrt(2 * p_demo * (1 - p_demo))
    
    G_demo_inc <- rbinom(N_pop, 2, p_demo)
    G_demo_dec <- rbinom(N_pop, 2, p_demo)
    
    G_raw[, demo_inc_idx] <- as.raw(G_demo_inc)
    G_raw[, demo_dec_idx] <- as.raw(G_demo_dec)
    
    demo_effect <- beta_demo * ((G_demo_inc - mu_d) / sd_d - (G_demo_dec - mu_d) / sd_d)
    gval <- gval + demo_effect
  } else {
    G_demo_inc <- NULL
    G_demo_dec <- NULL
  }
  
  # 6) add noise
  eps <- rnorm(N_pop, sd = sqrt(1 - h2))
  Y_pop <- as.numeric(scale(gval + eps))
  
  list(
    N_pop = N_pop,
    M = M,
    M_bg = M_bg,
    p = p,
    beta = beta,
    G_raw = G_raw,
    Y_pop = Y_pop,
    #demo handles
    demo_inc_idx = demo_inc_idx,
    demo_dec_idx = demo_dec_idx,
    h2 = h2,
    h2_bg = h2_bg,
    h2_demo = h2_demo
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
  r_mid <- (seq_len(B_eff) - 0.5) / B_eff # mid-quantile rank
  
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
# prep = prepped population phenotype
# N_samp = cohort size to draw
# q = threshold
# gamma = strength of skew (more weight on bins above q)
# eta = uniform mixing weight
# seed = RNG seed value
sample_cohort <- function(prep,
                          N_samp = 7000,
                          q = 0.50,
                          gamma = 5,
                          eta = 0.05,
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
  prob[keep_bins] <- eta * rep(1 / K, K) + (1 - eta) * p_keep
  
  cnt <- allocate_counts(N_samp, prob, prep$sizes)
  
  idx <- unlist(mapply(function(v, k) if (k > 0) sample(v, k, FALSE) else integer(0),
                       prep$idx_by_bin, cnt,
                       SIMPLIFY = FALSE, USE.NAMES = FALSE))
  y_samp <- prep$Y[idx]
  
  list(idx = idx, Y_samp = y_samp,
       mean = mean(y_samp), skew = skew3(y_samp),
       q = q, gamma = gamma, eta = eta)
}

### GWAS (linreg)
# pop = population phenotype object
# idx= indices of the cohort sampled
# chunk = which chunk to run
gwas_all_snps <- function(pop, idx, chunk = 200, ...) {
  stopifnot(!is.null(pop$G_raw), !is.null(pop$Y_pop))
  n <- length(idx)
  y <- pop$Y_pop[idx]
  y <- y - mean(y)
  Syy <- sum(y^2)
  
  M <- ncol(pop$G_raw)
  
  betahat <- rep(NA_real_, M)
  se      <- rep(NA_real_, M)
  maf     <- rep(NA_real_, M)
  mac     <- rep(NA_real_, M)
  af      <- rep(NA_real_, M)
  
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
    
    keep <- (mac_j >= 5) & (maf_j > 0)
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
  data.frame(
    snp = seq_len(M),
    af = af,
    maf = maf,
    mac = mac,
    betahat = betahat,
    se = se,
    z = z
  )
}

# ASH expected sign and mean sign-bias 
run_ash <- function(bhat, se) {
  if (!requireNamespace("ashr", quietly = TRUE)) {
    stop("Please install ashr: install.packages('ashr')")
  }
  ok <- is.finite(bhat) & is.finite(se) & (se > 0)
  bhat <- bhat[ok]
  se   <- se[ok]
  if (length(bhat) == 0L) return(numeric(0))
  
  s <- tryCatch({
    fit  <- ashr::ash(betahat = bhat, sebetahat = se,
                      method = "fdr",
                      mixcompdist = "normal",
                      optmethod = "mixSQP")
    lfsr <- ashr::get_lfsr(fit)
    sign(bhat) * (1 - 2 * lfsr)
  }, error = function(e) {
    z <- bhat / se
    sign(bhat) * (2 * pnorm(abs(z)) - 1)
  })
  
  s[is.finite(s)]
}

mean_sign_bias <- function(bhat, se) {
  s <- run_ash(bhat, se)
  if (length(s) == 0L) return(NA_real_)
  den <- sum(abs(s))
  if (!is.finite(den) || den <= .Machine$double.eps) return(NA_real_)
  sum(s) / den
}

.calc_p_from_z <- function(z) 2 * stats::pnorm(-abs(z))

.make_blocks <- function(M, block_size = 6) {
  B <- floor(M / block_size)
  split(seq_len(B * block_size), rep(seq_len(B), each = block_size))
}

sig_idx_from_gwas <- function(gwas,
                              select = c("block_min", "top_n"),
                              block_size = 6,
                              n_top = 1000,
                              z_col = "z") {
  select <- match.arg(select)
  z <- gwas[[z_col]]
  ok <- is.finite(z)
  p <- rep(Inf, length(z))
  p[ok] <- .calc_p_from_z(z[ok])
  
  if (select == "top_n") {
    ord <- order(p)
    ord <- ord[is.finite(p[ord])]
    return(head(ord, min(n_top, length(ord))))
  }
  
  blocks <- .make_blocks(nrow(gwas), block_size = block_size)
  idx <- integer(length(blocks))
  for (b in seq_along(blocks)) {
    block <- blocks[[b]]
    pb <- p[block]
    pb[!is.finite(pb)] <- Inf
    idx[b] <- block[which.min(pb)]
  }
  idx <- idx[is.finite(p[idx])]
  idx
}

.sample_block_minp <- function(p, blocks) {
  idx <- integer(length(blocks))
  for (b in seq_along(blocks)) {
    block <- blocks[[b]]
    pb <- p[block]
    pb[!is.finite(pb)] <- Inf
    idx[b] <- block[which.min(pb)]
  }
  idx
}

# Select SNPs (by block-min p or top_n) and compute ASH mean sign-bias
sign_bias_sig_ash <- function(gwas,
                              select = c("block_min", "top_n"),
                              block_size = 6,
                              n_top = 1000,
                              beta_col = "betahat",
                              se_col   = "se",
                              z_col    = "z") {
  select <- match.arg(select)
  
  req <- c(beta_col, se_col, z_col)
  miss <- setdiff(req, names(gwas))
  if (length(miss)) stop("gwas is missing required columns: ", paste(miss, collapse = ", "))
  
  bhat <- gwas[[beta_col]]
  se   <- gwas[[se_col]]
  z    <- gwas[[z_col]]
  
  ok <- is.finite(bhat) & is.finite(se) & (se > 0) & is.finite(z)
  if (!any(ok)) {
    return(data.frame(bias_sig = NA_real_, n_used = 0L))
  }
  
  p <- rep(Inf, length(z))
  p[ok] <- .calc_p_from_z(z[ok])
  
  if (select == "top_n") {
    ord <- order(p)
    ord <- ord[is.finite(p[ord])]
    idx <- head(ord, min(n_top, length(ord)))
  } else {
    if (nrow(gwas) < block_size) return(data.frame(bias_sig = NA_real_, n_used = 0L))
    blocks <- .make_blocks(nrow(gwas), block_size = block_size)
    idx <- .sample_block_minp(p, blocks)
    idx <- idx[is.finite(p[idx])]
  }
  
  data.frame(
    bias_sig = mean_sign_bias(bhat[idx], se[idx]),
    n_used   = length(idx)
  )
}

# Tune to hit a target skew
# grid search; returns best params and achieved skew
tune_params_for_skew <- function(prep,
                                 target_sk,
                                 N_samp,
                                 eta = 0.05,
                                 q_grid = c(0, 0.01, 0.05, 0.10,
                                            0.20, 0.30, 0.40,
                                            0.50, 0.60),
                                 gamma_grid = c(0, 1, 2, 3, 5, 8,
                                                13, 17, 20, 26, 35,
                                                50, 70, 100),
                                 seed = 1) {
  best <- NULL
  best_err <- Inf
  
  for (q in q_grid) {
    for (g in gamma_grid) {
      samp <- tryCatch(
        sample_cohort(prep,
                      N_samp = N_samp,
                      q = q,
                      gamma = g,
                      eta = eta,
                      seed = seed),
        error = function(e) NULL
      )
      if (is.null(samp)) next
      
      sk <- skew3(samp$Y_samp)
      err <- abs(sk - target_sk)
      
      if (is.finite(err) && err < best_err) {
        best_err <- err
        best <- list(q = q, gamma = g, skew = sk)
      }
    }
  }
  
  if (is.null(best)) {
    return(data.frame(target_sk = target_sk, q = NA_real_,
                      gamma = NA_real_,
                      skew_pilot = NA_real_, err = NA_real_))
  }
  
  data.frame(
    target_sk = target_sk,
    q = best$q,
    gamma = best$gamma,
    skew_pilot = best$skew,
    err = best_err
  )
}

#Get true sign bias
true_minor_sign <- function(af, beta_true) {
  s <- sign(beta_true)
  ifelse(af <= 0.5, s, -s)
}

#one cohort, sign bias
one_cohort_bias <- function(pop, prep,
                            target_sk, q, gamma,
                            N_samp,
                            eta = 0.05,
                            seed = 1,
                            gwas_chunk = 200,
                            mac_min = 0,
                            select = "block_min",
                            block_size = 6,
                            n_top = 1000,
                            give_sample = NULL) {
  
  samp <- if (!is.null(give_sample)) give_sample else
    sample_cohort(prep, N_samp = N_samp, q = q, gamma = gamma, eta = eta, seed = seed)
  
  gwas <- gwas_all_snps(pop, idx = samp$idx, chunk = gwas_chunk, mac_min = mac_min)
  
  idx_used <- sig_idx_from_gwas(gwas, select = select, block_size = block_size, n_top = n_top, z_col = "z")
  
  # estimated bias (ASH) on those SNPs
  sb <- sign_bias_sig_ash(gwas, select = select, block_size = block_size, n_top = n_top)
  
  # TRUE bias in this cohort, minor-allele oriented
  beta_true <- pop$beta[gwas$snp[idx_used]]
  s_true_minor <- true_minor_sign(gwas$af[idx_used], beta_true)
  true_bias_sample <- mean(s_true_minor)
  
  data.frame(
    target_sk = target_sk,
    skew_obs  = skew3(samp$Y_samp),
    q         = q,
    gamma     = gamma,
    bias_sig  = sb$bias_sig,
    true_bias_sample = true_bias_sample,
    n_used    = sb$n_used,
    seed      = seed
  )
}

#Skew vs sign bias
simulate_signbias_vs_skew <- function(pop, prep,
                                      targets = c(0, 0.5, 1, 1.5,
                                                  2, 2.5, 3, 3.5),
                                      n_rep = 20,
                                      N_samp = 10000,
                                      eta = 0.1,
                                      tune_seed = 1,
                                      q_grid = c(0, 0.01, 0.05, 0.10,
                                                 0.20, 0.30, 0.40, 0.50,
                                                 0.60),
                                      gamma_grid = c(0, 1, 2, 3, 5, 8, 13,
                                                     17, 20, 26, 35, 50, 70,
                                                     100),
                                      gwas_chunk = 200,
                                      mac_min = 0,
                                      select = "block_min",
                                      block_size = 6,
                                      n_top = 1000,
                                      ci_level = 0.95,
                                      base_seed = 1000) {
  
  # 1) Tune sampling params per target
  tuning <- do.call(rbind, lapply(seq_along(targets), function(i) {
    tune_params_for_skew(prep, target_sk = targets[i], N_samp = N_samp, eta = eta,
                         q_grid = q_grid, gamma_grid = gamma_grid, seed = tune_seed + i)
  }))
  
  # 2) Draw cohorts per target using tuned params
  reps <- do.call(rbind, lapply(seq_along(targets), function(i) {
    tsk <- targets[i]
    row <- tuning[tuning$target_sk == tsk, , drop = FALSE]
    
    if (!is.finite(row$q[1]) || !is.finite(row$gamma[1])) {
      return(data.frame(target_sk = tsk, skew_obs = NA_real_, q = NA_real_, gamma = NA_real_,
                        bias_sig = NA_real_, true_bias_sample = NA_real_, n_used = 0L, seed = NA_integer_))
    }
    
    do.call(rbind, lapply(seq_len(n_rep), function(r) {
      one_cohort_bias(pop, prep,
                      target_sk = tsk,
                      q = row$q[1],
                      gamma = row$gamma[1],
                      N_samp = N_samp,
                      eta = eta,
                      seed = base_seed + 10000L * i + r,
                      gwas_chunk = gwas_chunk,
                      mac_min = mac_min,
                      select = select,
                      block_size = block_size,
                      n_top = n_top)
    }))
  }))
  
  # 3) Summarize 
  alpha <- 1 - ci_level
  summary <- do.call(rbind, lapply(sort(unique(reps$target_sk)), function(tsk) {
    sub <- reps[reps$target_sk == tsk, , drop = FALSE]
    bs <- sub$bias_sig
    bs <- bs[is.finite(bs)]
    n <- length(bs)
    mean_bias <- if (n) mean(bs) else NA_real_
    sd_bias   <- if (n > 1) stats::sd(bs) else NA_real_
    se_mean   <- if (n > 1) sd_bias / sqrt(n) else NA_real_
    tcrit     <- if (n > 1) stats::qt(1 - alpha/2, df = n - 1) else NA_real_
    ci_low    <- mean_bias - tcrit * se_mean
    ci_high   <- mean_bias + tcrit * se_mean
    data.frame(target_sk = tsk, n_rep = n,
               mean_skew_obs = mean(sub$skew_obs, na.rm = TRUE),
               mean_bias = mean_bias, sd_bias = sd_bias,
               se_mean = se_mean, ci_low = ci_low, ci_high = ci_high)
  }))
  
  list(tuning = tuning, reps = reps, summary = summary)
}

## Common params
N_pop     <- .as_int(.get_value(args, "--N_pop", "2000000"), "--N_pop")
pop_seed  <- .as_int(.get_value(args, "--pop_seed", "1"), "--pop_seed")
chunk     <- .as_int(.get_value(args, "--chunk", "200"), "--chunk")
B_bins    <- .as_int(.get_value(args, "--B_bins", "200"), "--B_bins")

out_sample <- .get_value(args, "--out_sample", "sample.rds")
out_sim    <- .get_value(args, "--out_sim", "out_reps.txt")

## Build population + prep once
cat("Building population...\n")
pop <- build_population(N_pop = N_pop, seed = pop_seed, chunk = chunk)
cat("Prepping bins...\n")
prep <- prep_bins(pop$Y_pop, B = B_bins)

## --single: write sample.rds for single cohort
if (run_single) {
  N_samp      <- .as_int(.get_value(args, "--N_samp", "10000"), "--N_samp")
  q           <- .as_num(.get_value(args, "--q", "0"), "--q")
  gamma       <- .as_num(.get_value(args, "--gamma", "19.2"), "--gamma")
  eta         <- .as_num(.get_value(args, "--eta", "0.1"), "--eta")
  sample_seed <- .as_int(.get_value(args, "--sample_seed", "199"), "--sample_seed")
  panels      <- .as_bool(.get_value(args, "--panels", "false"), "--panels")
  
  ## cohort
  demo_cohort_args <- sample_cohort(prep,
                                    N_samp = N_samp,
                                    q = q,
                                    gamma = gamma,
                                    eta = eta,
                                    seed = sample_seed)
  
  Y_pop  <- pop$Y_pop
  Y_demo <- demo_cohort_args$Y_samp
  
  ## Demo SNP genotypes from the user-specified cohort
  if (is.na(pop$demo_inc_idx) || is.na(pop$demo_dec_idx)) {
    stop("Population was built with add_demo=FALSE; cannot compute G_inc/G_dec needed for plotting.")
  }
  
  G_inc_raw <- as.integer(pop$G_raw[demo_cohort_args$idx, pop$demo_inc_idx, drop = TRUE])
  G_dec_raw <- as.integer(pop$G_raw[demo_cohort_args$idx, pop$demo_dec_idx, drop = TRUE])
  
  ## Orient inc/dec by fitted slopes (exactly like your original Panel B code)
  b_inc <- stats::coef(stats::lm(Y_demo ~ G_inc_raw))[2]
  b_dec <- stats::coef(stats::lm(Y_demo ~ G_dec_raw))[2]
  if (b_inc >= b_dec) {
    G_inc <- G_inc_raw
    G_dec <- G_dec_raw
  } else {
    G_inc <- G_dec_raw
    G_dec <- G_inc_raw
  }
  
  ## gwas_cohort, panels
  if (isTRUE(panels)) {
    demo_cohort_panel <- sample_cohort(prep,
                                       N_samp = 10000,
                                       q = 0.5,
                                       gamma = 3.39,
                                       eta = 0.1,
                                       seed = 1)
    gwas_cohort <- gwas_all_snps(pop, idx = demo_cohort_panel$idx, chunk = chunk)
  } else {
    gwas_cohort <- gwas_all_snps(pop, idx = demo_cohort_args$idx, chunk = chunk)
  }
  
  sims_run_f3 <- list(Y_pop, Y_demo, G_inc, G_dec, gwas_cohort, pb = pop$beta)
  saveRDS(sims_run_f3, file = out_sample)
  cat("Wrote sample object to: ", out_sample, "\n", sep = "")
}

## --sim: write out_reps as plain text to --out_sim
if (run_sim) {
  n_rep <- .as_int(.get_value(args, "--n_rep", "20"), "--n_rep")
  
  targets <- .parse_num_list(.get_value(args, "--targets", ""))
  if (is.null(targets)) {
    t0 <- .as_num(.get_value(args, "--targets_start", "0"), "--targets_start")
    t1 <- .as_num(.get_value(args, "--targets_end", "4"), "--targets_end")
    ts <- .as_num(.get_value(args, "--targets_step", "0.1"), "--targets_step")
    targets <- seq(t0, t1, by = ts)
  }

  out <- simulate_signbias_vs_skew(pop, prep,
                                   targets = targets,
                                   n_rep = n_rep,
                                   N_samp = 10000,
                                   eta = 0.1)
  
  out_reps <- out$reps
  
  ## Write
  write.table(out_reps, file = out_sim, sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  cat("Wrote out_reps table to: ", out_sim, "\n", sep = "")
}

cat("Done.\n")

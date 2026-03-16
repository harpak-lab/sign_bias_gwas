############################################################
## Simulation of skewed cohorts, scheme B
##
## A quantitative trait is simulated in a population with three
## latent phenotype modes (left, center, right), with mode- and
## sign-dependent allele frequencies. Cohorts are then sampled
## using a two-bump scheme with fixed center and right kernels,
## while only the right-mode sampling fraction varies. 
##
## Parameters:
##   N_POP         population size
##   N_VARIANTS    number of simulated SNPs
##   H2            target SNP heritability
##   N_COHORT      cohort sample size
##   TAIL_FRAC     fraction in each population tail mode
##   MAF_SHIFT     mode-dependent allele-frequency shift
##   RARE_MAF_MAX  MAF threshold for evaluation
##   RARE_MAC_MIN  MAC threshold for evaluation
##   BLOCK_SIZE    pseudo-block size for min-p winner selection
##   N_REPS        replicates per target skew condition
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ashr)
library(cowplot)

#Parameters, input
N_POP         <- 2000000
N_VARIANTS    <- 4000
H2            <- 0.6
MAF_MIN_POP   <- 1e-4
MAF_MAX_POP   <- 0.1
POP_SEED      <- 42
STORE_G       <- TRUE
TAIL_FRAC     <- 0.1
MAF_SHIFT     <- 0.5

RARE_MAF_MAX  <- 0.01
RARE_MAC_MIN  <- 1
BLOCK_SIZE    <- 6

N_COHORT      <- 10000
N_REPS        <- 20
SIM_SEED      <- 1
VERBOSE       <- TRUE
ADD_SKEW0     <- FALSE
N_REP0        <- 20

CENTER_LAMBDA <- 200
RIGHT_LAMBDA  <- 2
BAND_SD_CENTER <- 0.1
BAND_SD_RIGHT  <- 0.037
WEIGHT_ON      <- "y"

RIGHT_FRAC_MAP <- c(
  `0`  = 0.00,
  `1`  = 0.276,
  `2`  = 0.146,
  `3`  = 0.0835,
  `4`  = 0.0523,
  `5`  = 0.0355,
  `6`  = 0.0254,
  `7`  = 0.0189,
  `8`  = 0.0146,
  `9`  = 0.01165,
  `10` = 0.0094
)

RARE_MAF_MAX <- RARE_MAF_MAX
RARE_MAC_MIN <- RARE_MAC_MIN
BLOCK_SIZE   <- BLOCK_SIZE

log_msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) cat(sprintf(...), "\n")
}

skew3 <- function(x){
  m <- mean(x)
  s <- sd(x)
  if (s == 0) 0 else mean(((x - m) / s)^3)
}

lowest_p_per6_inorder <- function(df, block_size = BLOCK_SIZE){
  if (nrow(df) == 0) return(df[0, , drop = FALSE])
  df <- df %>% arrange(variant)
  block_id <- ((df$variant - 1L) %/% as.integer(block_size)) + 1L
  df %>%
    mutate(block_id = block_id) %>%
    group_by(block_id) %>%
    slice_min(order_by = p, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-block_id)
}

# 1) Population builder (3-mode; spread-out via genotype comp)
build_population_gmix_psbased <- function(
    n_pop      = 300000,
    n_variants = 2000,      # must be even
    h2 = 0.6,      # equal magnitude across all SNPs
    maf_min    = 1e-4,
    maf_max    = 0.1,
    seed       = 42,
    store_G    = TRUE,
    tail_frac  = 0.05,
    maf_shift  = 0.7
){
  stopifnot(n_variants %% 2 == 0)
  stopifnot(tail_frac > 0 && tail_frac < 0.5)
  stopifnot(is.finite(maf_shift) && maf_shift >= 0 && maf_shift < 1)
  
  set.seed(seed)
  
  M2 <- n_variants / 2
  
  maf_half <- 10^(runif(M2, log10(maf_min), log10(maf_max)))
  mafs <- rep(maf_half, times = 2)
  
  sum_varG_bg <- sum(2 * mafs * (1 - mafs))
  effect_size <- sqrt(h2 / sum_varG_bg)
  
  effects <- c(rep(+effect_size, M2), rep(-effect_size, M2))
  sign_true <- sign(effects)
  pair_id <- rep(seq_len(M2), times = 2)
  
  nL <- as.integer(round(n_pop * tail_frac))
  nR <- as.integer(round(n_pop * tail_frac))
  nC <- n_pop - nL - nR
  if (nC <= 0) stop("tail_frac too large: center mode would be empty.")
  
  z <- rep("center", n_pop)
  z[sample.int(n_pop, nL)] <- "left"
  remaining <- which(z == "center")
  z[sample(remaining, nR)] <- "right"
  z <- factor(z, levels = c("left","center","right"))
  
  G_raw <- if (store_G) matrix(as.raw(0), nrow = n_pop, ncol = n_variants) else NULL
  g <- numeric(n_pop)
  
  for (j in seq_len(n_variants)) {
    p0 <- mafs[j]
    
    if (effects[j] > 0) {
      pL <- p0 * (1 - maf_shift)
      pC <- p0
      pR <- p0 * (1 + maf_shift)
    } else {
      pL <- p0 * (1 + maf_shift)
      pC <- p0
      pR <- p0 * (1 - maf_shift)
    }
    
    pL <- min(max(pL, 1e-8), 1 - 1e-8)
    pC <- min(max(pC, 1e-8), 1 - 1e-8)
    pR <- min(max(pR, 1e-8), 1 - 1e-8)
    
    pvec <- ifelse(z == "left", pL, ifelse(z == "right", pR, pC))
    gj <- as.integer(rbinom(n_pop, size = 2, prob = pvec))
    
    if (store_G) G_raw[, j] <- as.raw(gj)
    g <- g + effects[j] * gj
  }
  
  y <- g + rnorm(n_pop, mean = 0, sd = 1-h2)
  
  snp_meta <- data.frame(
    snp_id    = seq_len(n_variants),
    pair_id   = pair_id,
    maf       = mafs,
    effect    = effects,
    sign_true = sign_true
  )
  
  mode_sizes <- as.integer(table(z)); names(mode_sizes) <- levels(z)
  
  cat(
    "Population y: mean=", mean(y),
    " sd=", sd(y),
    " skew=", skew3(y), "\n",
    "Mode sizes:", paste(levels(z), mode_sizes, collapse="  "), "\n",
    "mean(sign(true_effect)) across SNPs =", mean(sign(effects)), "\n",
    "cor(g, y)=", cor(g, y), "\n",
    "Mode means of g:", paste(levels(z), round(tapply(g, z, mean), 4), collapse="  "), "\n",
    "Mode means of y:", paste(levels(z), round(tapply(y, z, mean), 4), collapse="  "), "\n"
  )
  
  list(
    y = y,
    g = g,
    z = z,
    G_raw = G_raw,
    snp_meta = snp_meta,
    tail_frac = tail_frac,
    maf_shift = maf_shift
  )
}

# 2) Cohort sampler: 2-bump (tight center + tight right bump)
sample_cohort_2bump <- function(pop, n,
                                right_frac = 0.05,
                                center_lambda = 200,
                                right_lambda  = 2,
                                band_sd_center = 0.1,
                                band_sd_right  = 0.03,
                                seed = 1,
                                weight_on = c("y","g")) {
  
  set.seed(seed)
  weight_on <- match.arg(weight_on)
  
  x <- switch(weight_on,
              y = pop$y,
              g = pop$g)
  
  idxC <- which(pop$z == "center")
  idxR <- which(pop$z == "right")
  
  nR <- round(n * right_frac)
  nC <- n - nR
  
  if (nC > length(idxC)) stop("Not enough center-mode individuals.")
  if (nR > length(idxR)) stop("Not enough right-mode individuals.")
  
  xC <- x[idxC]
  muC <- median(xC)
  sdC <- sd(xC) + 1e-12
  keepC <- abs(xC - muC) <= band_sd_center * sdC
  idxC2 <- idxC[keepC]
  xC2 <- x[idxC2]
  if (nC > length(idxC2)) stop("Center band too tight; increase band_sd_center or lower n.")
  wC <- exp(-center_lambda * ((xC2 - muC) / sdC)^2)
  
  xR <- x[idxR]
  muR <- median(xR)
  sdR <- sd(xR) + 1e-12
  keepR <- abs(xR - muR) <= band_sd_right * sdR
  idxR2 <- idxR[keepR]
  xR2 <- x[idxR2]
  if (nR > length(idxR2)) stop("Right band too tight; increase band_sd_right or lower n.")
  wR <- exp(-right_lambda * ((xR2 - muR) / sdR)^2)
  
  pickC <- sample(idxC2, nC, replace = FALSE, prob = wC)
  pickR <- if (nR > 0) sample(idxR2, nR, replace = FALSE, prob = wR) else integer(0)
  
  idx <- sample(c(pickC, pickR))
  y_s <- pop$y[idx]
  
  list(
    idx = idx,
    y = y_s,
    mean = mean(y_s),
    sd = sd(y_s),
    skew = skew3(y_s),
    right_frac = right_frac,
    weight_on = weight_on
  )
}

# 3) GWAS: fast single-SNP OLS
run_gwas_ols_fast <- function(y, G_int, effects_true){
  stopifnot(length(y) == nrow(G_int), length(effects_true) == ncol(G_int))
  n <- length(y)
  m <- ncol(G_int)
  
  y0 <- as.numeric(y) - mean(y0 <- as.numeric(y))
  Syy <- sum(y0^2)
  
  ac  <- colSums(G_int)
  af  <- ac / (2 * n)
  maf <- pmin(af, 1 - af)
  mac <- pmin(ac, 2*n - ac)
  
  gbar <- ac / n
  G2sum <- colSums(G_int * G_int)
  Sxx <- G2sum - n * (gbar^2)
  
  Sxy <- as.numeric(crossprod(G_int, y0))
  
  ok <- is.finite(Sxx) & (Sxx > 0)
  
  beta <- rep(NA, m)
  beta[ok] <- Sxy[ok] / Sxx[ok]
  
  SSE <- rep(NA, m)
  SSE[ok] <- pmax(Syy - beta[ok] * Sxy[ok], 0)
  
  sigma2_hat <- rep(NA, m)
  sigma2_hat[ok] <- SSE[ok] / (n - 2)
  
  se <- rep(NA, m)
  se[ok] <- sqrt(sigma2_hat[ok] / Sxx[ok])
  
  tstat <- rep(NA, m)
  pval  <- rep(NA, m)
  tstat[ok] <- beta[ok] / se[ok]
  pval[ok]  <- 2 * pt(-abs(tstat[ok]), df = n - 2)
  
  tibble(
    variant = seq_len(m),
    true_effect = effects_true,
    effect_direction = factor(ifelse(effects_true > 0, "increasing", "decreasing"),
                              levels = c("increasing","decreasing")),
    af_sample  = af,
    maf_sample = maf,
    ac_sample  = ac,
    mac_sample = mac,
    beta = beta,
    se = se,
    t = tstat,
    p = pval,
    n = n
  )
}

# 4) ASH per-SNP sign bias: s = Pr(beta>0) - Pr(beta<0)
ash_signbias_per_snp <- function(beta_hat, se_hat,
                                 mixcompdist = "normal",
                                 method = "fdr") {
  n <- length(beta_hat)
  sb <- rep(NA, n)
  
  ok <- is.finite(beta_hat) & is.finite(se_hat) & se_hat > 0
  if (!any(ok)) return(sb)
  
  sb_ok <- tryCatch({
    fit <- ashr::ash(betahat = beta_hat[ok],
                     sebetahat = se_hat[ok],
                     mixcompdist = mixcompdist,
                     method = method)
    
    if (!is.null(fit$result) &&
        all(c("PositiveProb","NegativeProb") %in% colnames(fit$result))) {
      as.numeric(fit$result$PositiveProb - fit$result$NegativeProb)
    } else {
      p_ge0 <- as.numeric(ashr::get_posterior_prob(fit, l = 0,    u = Inf))  # Pr(beta>=0)
      p_le0 <- as.numeric(ashr::get_posterior_prob(fit, l = -Inf, u = 0))    # Pr(beta<=0)
      as.numeric(p_ge0 - p_le0)                                             # == Pr(beta>0)-Pr(beta<0)
    }
  }, error = function(e){
    z <- beta_hat[ok] / se_hat[ok]
    sign(z) * (2 * pnorm(abs(z)) - 1)
  })
  
  sb[ok] <- sb_ok
  sb
}

# 5) One cohort
one_cohort_metrics <- function(pop, coh, maf_max = RARE_MAF_MAX, mac_min = RARE_MAC_MIN){
  idx <- coh$idx
  y   <- coh$y
  
  G_int <- matrix(as.integer(pop$G_raw[idx, , drop = FALSE]),
                  nrow = length(idx), ncol = ncol(pop$G_raw))
  
  effects_true <- pop$snp_meta$effect
  gwas <- run_gwas_ols_fast(y, G_int, effects_true)
  
  # ASH per SNP (all SNPs)
  gwas$ash_s <- ash_signbias_per_snp(gwas$beta, gwas$se)
  
  keep <- is.finite(gwas$maf_sample) & (gwas$maf_sample > 0) & (gwas$maf_sample <= maf_max) &
    is.finite(gwas$mac_sample) & (gwas$mac_sample >= mac_min) &
    is.finite(gwas$true_effect) & (gwas$true_effect != 0) &
    is.finite(gwas$beta) & is.finite(gwas$se) & (gwas$se > 0) &
    is.finite(gwas$p)
  
  # TRUE 
  true_bias <- if (any(keep)) mean(sign(gwas$true_effect[keep])) else NA
  sub <- gwas %>% filter(keep & is.finite(ash_s))
  win <- lowest_p_per6_inorder(sub, block_size = BLOCK_SIZE)
  
  if (nrow(win) == 0) {
    ash_bias <- NA
    n_selected <- 0
  } else {
    s <- win$ash_s
    num <- sum(s, na.rm = TRUE)
    den <- sum(abs(s), na.rm = TRUE)
    ash_bias <- if (is.finite(den) && den > 0) num / den else NA
    n_selected <- nrow(win)
  }
  # -------------------------------------------------------------------------
  
  # ratios (inc/dec) 
  if (!any(keep)) {
    ratio_absb <- ratio_se <- ratio_t2 <- NA
  } else {
    gg <- gwas[keep, ]
    inc <- gg$true_effect > 0
    dec <- gg$true_effect < 0
    if (!any(inc) || !any(dec)) {
      ratio_absb <- ratio_se <- ratio_t2 <- NA
    } else {
      ratio_absb <- mean(abs(gg$beta[inc]), na.rm=TRUE) / mean(abs(gg$beta[dec]), na.rm=TRUE)
      ratio_se   <- mean(gg$se[inc], na.rm=TRUE)       / mean(gg$se[dec], na.rm=TRUE)
      ratio_t2   <- mean((gg$t[inc])^2, na.rm=TRUE)    / mean((gg$t[dec])^2, na.rm=TRUE)
    }
  }
  
  tibble(
    realized_skew = coh$skew,
    true_bias = true_bias,
    ash_bias  = ash_bias,
    n_selected = n_selected,
    ratio_absb = ratio_absb,
    ratio_se   = ratio_se,
    ratio_t2   = ratio_t2
  )
}


# 6) Experiment run
run_experiment <- function(
    pop,
    n_cohort = 10000,
    n_reps = 20,
    right_frac_map = c(
      `1`  = 0.28,
      `2`  = 0.14,
      `3`  = 0.08,
      `4`  = 0.05,
      `5`  = 0.035,
      `6`  = 0.025,
      `7`  = 0.0185,
      `8`  = 0.0137,
      `9`  = 0.0105,
      `10` = 0.0084
    ),
    sampler_fixed = list(
      center_lambda  = 200,
      right_lambda   = 2,
      band_sd_center = 0.1,
      band_sd_right  = 0.03,
      weight_on      = "y"
    ),
    maf_max = RARE_MAF_MAX,
    mac_min = RARE_MAC_MIN,
    sim_seed = 10001,
    verbose = TRUE,
    add_skew0 = TRUE,
    n_rep0 = 20
){
  stopifnot(!is.null(pop$G_raw), !is.null(pop$y), !is.null(pop$snp_meta))
  
  skew_targets <- sort(as.integer(names(right_frac_map)))
  rows <- list(); ir <- 0
  
  log_msg("=== START ===", verbose=verbose)
  log_msg("Cohort n=%d | SNPs=%d | reps=%d | targets=%s",
          n_cohort, ncol(pop$G_raw), n_reps, paste(skew_targets, collapse=","), verbose=verbose)
  log_msg("MAF threshold = %.3g | MAC min = %d", maf_max, mac_min, verbose=verbose)
  
  if (isTRUE(add_skew0)) {
    log_msg("\n[BASELINE skew=0] n_rep0=%d", n_rep0, verbose=verbose)
    
    for (r in seq_len(n_rep0)) {
      this_seed <- sim_seed + 900000 + r
      
      coh0 <- do.call(sample_cohort_2bump, c(
        list(pop = pop, n = n_cohort, right_frac = 0, seed = this_seed),
        modifyList(sampler_fixed, list(center_lambda = 0, band_sd_center = 10))
      ))
      
      met <- one_cohort_metrics(pop, coh0, maf_max = maf_max, mac_min = mac_min)
      
      ir <- ir + 1L
      rows[[ir]] <- tibble(
        target_skew = 0,
        rep = r,
        right_frac = 0,
        realized_skew = met$realized_skew,
        true_bias = met$true_bias,
        ash_bias  = met$ash_bias,
        n_selected = met$n_selected,
        ratio_absb = met$ratio_absb,
        ratio_se   = met$ratio_se,
        ratio_t2   = met$ratio_t2
      )
    }
  }
  
  for (tg in skew_targets){
    rf <- unname(right_frac_map[as.character(tg)])
    log_msg("\n[TARGET %d] right_frac=%.5f", tg, rf, verbose=verbose)
    
    for (r in seq_len(n_reps)){
      this_seed <- sim_seed + 100000*tg + r
      
      coh <- do.call(sample_cohort_2bump, c(
        list(pop = pop, n = n_cohort, right_frac = rf, seed = this_seed),
        sampler_fixed
      ))
      
      met <- one_cohort_metrics(pop, coh, maf_max = maf_max, mac_min = mac_min)
      
      ir <- ir + 1L
      rows[[ir]] <- tibble(
        target_skew = as.integer(tg),
        rep = r,
        right_frac = rf,
        realized_skew = met$realized_skew,
        true_bias = met$true_bias,
        ash_bias  = met$ash_bias,
        n_selected = met$n_selected,
        ratio_absb = met$ratio_absb,
        ratio_se   = met$ratio_se,
        ratio_t2   = met$ratio_t2
      )
    }
  }
  
  df <- bind_rows(rows)
  
  log_msg("\n=== END ===", verbose=verbose)
  log_msg("rows=%d", nrow(df), verbose=verbose)
  
  df
}

# 7) Summaries + plots
summarize_bias_only <- function(df, ci_level = 0.95, drop_zero = TRUE){
  alpha <- 1 - ci_level
  
  if (drop_zero) df <- df %>% filter(target_skew != 0)
  
  df %>%
    group_by(target_skew) %>%
    summarise(
      n_rep = sum(is.finite(true_bias) & is.finite(ash_bias)),
      
      mean_true = mean(true_bias, na.rm = TRUE),
      sd_true   = sd(true_bias, na.rm = TRUE),
      se_true   = sd_true / sqrt(sum(is.finite(true_bias))),
      
      mean_ash = mean(ash_bias, na.rm = TRUE),
      sd_ash   = sd(ash_bias, na.rm = TRUE),
      se_ash   = sd_ash / sqrt(sum(is.finite(ash_bias))),
      
      tcrit = ifelse(n_rep > 1, qt(1 - alpha/2, df = n_rep - 1), NA),
      
      lo_true = mean_true - tcrit * se_true,
      hi_true = mean_true + tcrit * se_true,
      
      lo_ash  = mean_ash  - tcrit * se_ash,
      hi_ash  = mean_ash  + tcrit * se_ash,
      
      .groups = "drop"
    )
}

plot_bias_panelE_style <- function(bsum){
  pdat <- bind_rows(
    bsum %>% transmute(x = target_skew,
                       series = "Cohort true sign bias",
                       mean = mean_true, ymin = lo_true, ymax = hi_true, shape = "est"),
    bsum %>% transmute(x = target_skew,
                       series = "Estimated sign bias",
                       mean = mean_ash,  ymin = lo_ash,  ymax = hi_ash,  shape = "est"),
    bsum %>% transmute(x = target_skew,
                       series = "Population true sign bias",
                       mean = 0,  ymin = 0,  ymax = 0,  shape = "est")
  )
  
  ggplot(pdat, aes(x = x, y = mean, color = series, group = series)) +
    geom_line(linewidth = 0.7) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.5) +
    geom_point(aes(shape = shape), size = 2.0) +
    scale_color_manual(values = c(
      "Cohort true sign bias" = "purple4",
      "Estimated sign bias"   = "black",
      "Population true sign bias"   = "#D5B60A"
    )) +
    scale_x_continuous(breaks = pretty(bsum$target_skew),
                       limits = range(bsum$target_skew)) +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_text(size = 12),
      legend.position = "none"
    )
}

summarize_ratios <- function(df, ci_level = 0.95){
  alpha <- 1 - ci_level
  df %>%
    group_by(target_skew) %>%
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
}

plot_ratios_original_style <- function(diag_summ){
  diag_plot_df <- dplyr::bind_rows(
    diag_summ %>% dplyr::transmute(target_skew, metric = "|beta| mean ratio (inc/dec)",
                                   mean = mean_absb, lo = lo_absb, hi = hi_absb),
    diag_summ %>% dplyr::transmute(target_skew, metric = "SE mean ratio (inc/dec)",
                                   mean = mean_se, lo = lo_se, hi = hi_se),
    diag_summ %>% dplyr::transmute(target_skew, metric = "T^2 mean ratio (inc/dec)",
                                   mean = mean_t2, lo = lo_t2, hi = hi_t2)
  ) %>%
    dplyr::mutate(
      target_sk_f = factor(target_skew, levels = rev(sort(unique(target_skew))))
    )
  
  ggplot(diag_plot_df, aes(x = mean, y = target_sk_f)) +
    geom_vline(xintercept = 1, linewidth = 0.8, color = "#D5B60A") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.6) +
    geom_point(size = 2.2) +
    facet_wrap(~ metric, ncol = 1, scales = "free_x") +
    theme_classic(base_size = 12) +
    labs(
      x = "Ratio (increasing / decreasing)",
      y = "Target skew",
      title = paste0("GWAS diagnostic ratios (MAF \u2264 ", format(RARE_MAF_MAX, scientific=TRUE), ")")
    ) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
}

# 8) RUN

pop1 <- build_population_gmix_psbased(
  n_pop      = N_POP,
  n_variants = N_VARIANTS,
  maf_min    = MAF_MIN_POP,
  maf_max    = MAF_MAX_POP,
  seed       = POP_SEED,
  store_G    = STORE_G,
  tail_frac  = TAIL_FRAC,
  maf_shift  = MAF_SHIFT
)

right_frac_map <- RIGHT_FRAC_MAP

sampler_fixed <- list(
  center_lambda  = CENTER_LAMBDA,
  right_lambda   = RIGHT_LAMBDA,
  band_sd_center = BAND_SD_CENTER,
  band_sd_right  = BAND_SD_RIGHT,
  weight_on      = WEIGHT_ON
)

res_df <- run_experiment(
  pop1,
  n_cohort = N_COHORT,
  n_reps = N_REPS,
  right_frac_map = right_frac_map,
  sampler_fixed = sampler_fixed,
  maf_max = RARE_MAF_MAX,
  mac_min = RARE_MAC_MIN,
  sim_seed = SIM_SEED,
  verbose = VERBOSE,
  add_skew0 = ADD_SKEW0,
  n_rep0 = N_REP0
)

dens_pop <- density(pop1$y)
dd_pop <- data.frame(x=dens_pop$x,
                     y=dens_pop$y)

mode_centers <- data.frame(
  z  = levels(pop1$z),
  x0 = as.numeric(tapply(pop1$y, pop1$z, mean))
)

pD <- ggplot(dd_pop, aes(x,y)) +
  geom_area(fill="#FDE992", color=NA) +
  geom_line(color="black", linewidth=0.6) +
  geom_vline(
    data = mode_centers,
    aes(xintercept = x0),
    linewidth = 0.8,
    linetype = 2
  ) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10))

coh <- sample_cohort_2bump(pop1, n=N_COHORT,
                           center_lambda  = CENTER_LAMBDA,
                           right_lambda   = RIGHT_LAMBDA,
                           band_sd_center = BAND_SD_CENTER,
                           band_sd_right  = BAND_SD_RIGHT,
                           weight_on      = WEIGHT_ON,
                           right_frac =  0.0094)

Y_demo <- coh$y

dens_samp <- density(Y_demo, adjust=3)
dd_samp <- data.frame(x=dens_samp$x,
                      y=dens_samp$y)

pE <- ggplot(dd_samp, aes(x,y)) +
  geom_area(fill="#EFEFEF", color=NA) +
  geom_line(color="black", linewidth=0.6) +
  theme_classic() +
  theme(axis.title.y     = element_blank(),
        axis.title.x     = element_blank(),
        axis.text.x      = element_text(size = 10),
        axis.text.y      = element_text(size = 10))

bias_sum <- summarize_bias_only(res_df, drop_zero = FALSE)
print(bias_sum)

pF <- plot_bias_panelE_style(bias_sum)
pF <- pF + ylim(-0.1,0.4) + xlim(0,10.1)

diag_summ <- summarize_ratios(res_df)
print(diag_summ)
pDiag <- plot_ratios_original_style(diag_summ)
print(pDiag)

res_df %>%
  group_by(target_skew) %>%
  summarize(mean_realized_skew = mean(realized_skew, na.rm=TRUE),
            sd_realized_skew   = sd(realized_skew, na.rm=TRUE),
            .groups="drop") %>%
  print()

pdat <- bind_rows(
  bias_sum %>% transmute(
    x = target_skew,
    series = "Cohort true sign bias",
    mean = mean_true, ymin = lo_true, ymax = hi_true
  ),
  bias_sum %>% transmute(
    x = target_skew,
    series = "Estimated sign bias",
    mean = mean_ash,  ymin = lo_ash,  ymax = hi_ash
  ),
  bias_sum %>% transmute(
    x = target_skew,
    series = "Population true sign bias",
    mean = 0, ymin = NA, ymax = NA
  )
)

p_y <- dplyr::filter(pdat, series == "Population true sign bias")
p_p <- dplyr::filter(pdat, series == "Cohort true sign bias")
p_b <- dplyr::filter(pdat, series == "Estimated sign bias")


pF <- ggplot() +
  ## 1) Yellow (bottom)
  geom_line(data = p_y, aes(x = x, y = mean, group = 1),
            linewidth = 0.7, color = "#D5B60A") +
  geom_point(data = p_y, aes(x = x, y = mean),
             shape = 16, size = 2.0, color = "#D5B60A") +
  
  ## 2) Purple (middle, ABOVE yellow)
  geom_line(data = p_p, aes(x = x, y = mean, group = 1),
            linewidth = 0.7, color = "purple4") +
  geom_errorbar(
    data = dplyr::filter(p_p, !is.na(ymin) & !is.na(ymax)),
    aes(x = x, ymin = ymin, ymax = ymax),
    width = 0.15, linewidth = 0.5, color = "purple4"
  ) +
  geom_point(data = p_p, aes(x = x, y = mean),
             shape = 16, size = 2.0, color = "purple4") +
  
  ## 3) Black (top)
  geom_line(data = p_b, aes(x = x, y = mean, group = 1),
            linewidth = 0.7, color = "black") +
  geom_errorbar(
    data = dplyr::filter(p_b, !is.na(ymin) & !is.na(ymax)),
    aes(x = x, ymin = ymin, ymax = ymax),
    width = 0.15, linewidth = 0.5, color = "black"
  ) +
  geom_point(data = p_b, aes(x = x, y = mean),
             shape = 16, size = 2.0, color = "black") +
  
  scale_x_continuous(
    breaks = pretty(bias_sum$target_skew),
    limits = range(bias_sum$target_skew)
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_text(size = 10),
    legend.position = "none"
  ) +
  xlim(1, 10.1) +
  ylim(-0.05,0.4)

#### plot

bottom_fig3<-plot_grid(pD,pE,pF, nrow=1)
bottom_fig3
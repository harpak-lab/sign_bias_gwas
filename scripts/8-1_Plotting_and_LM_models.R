##################################################
## Plotting and summary analyses
##################################################
setwd("~/updated_results/")

#libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(grid)
library(gridExtra)
library(patchwork)
library(scales)
library(forcats)

##################################################
## Main text 

########## FIG 1

### load data
h2_ME_all <- fread("UKB_vs_AOU_h2_Me.csv", sep=",")
gcorr <- fread("gcorr.csv")

### prepare data
h2_ME <- h2_ME_all
h2_ME <- h2_ME %>%
  pivot_wider(
    names_from  = database,
    values_from = c(h2, h2_se, Me, Me_se, log10_Me, log10_Me_SE),
    names_sep   = "_"
  )

### calculate ratios and CIs
h2_ME <- h2_ME %>%
  mutate(
    ratio_h2 = h2_AoU / h2_UKB,
    ratio_ME = Me_AoU / Me_UKB,
    
    se_log_ratio_h2 = sqrt( (h2_se_AoU / h2_AoU)^2 + (h2_se_UKB / h2_UKB)^2 ),
    se_log_ratio_ME = sqrt( (Me_se_AoU / Me_AoU)^2 + (Me_se_UKB / Me_UKB)^2 ),
    
    min_h2 = exp(log(ratio_h2) - 1.96 * se_log_ratio_h2),
    max_h2 = exp(log(ratio_h2) + 1.96 * se_log_ratio_h2),
    min_ME = exp(log(ratio_ME) - 1.96 * se_log_ratio_ME),
    max_ME = exp(log(ratio_ME) + 1.96 * se_log_ratio_ME),
    
    se_ratio_h2 = ratio_h2 * se_log_ratio_h2,   
    min_h2_lin  = ratio_h2 - 1.96 * se_ratio_h2,
    max_h2_lin  = ratio_h2 + 1.96 * se_ratio_h2
  )

h2_ME<-merge(h2_ME, gcorr, by="Trait")

h2_ME$min_rg <- h2_ME$rg - 1.96 * h2_ME$SE
h2_ME$max_rg <- h2_ME$rg + 1.96 * h2_ME$SE

### Make plot
levs <- with(h2_ME, levels(reorder(factor(Trait), -ratio_h2)))
h2_ME$Trait_ord <- factor(h2_ME$Trait, levels = levs)

## add stripe
n <- length(levs)
stripe_df <- data.frame(
  y    = seq_len(n),
  ymin = seq_len(n) - 0.5,
  ymax = seq_len(n) + 0.5
)
stripe_df <- subset(stripe_df, y %% 2 == 0)

## Make panel function
make_panel <- function(dat, xvar, xmin, xmax,
                       show_y     = TRUE,
                       point_col  = "#B56727",
                       log_scale  = FALSE,
                       log_base   = 10,     # <-- base-2 fits ratio intuition (…0.5,1,2,4)
                       tick_n     = 5,     # <-- target number of major ticks
                       ref_at     = 1) {
  
  dat <- dat |>
    dplyr::mutate(.x = .data[[xvar]], .xmin = .data[[xmin]], .xmax = .data[[xmax]])
  
  if (isTRUE(log_scale)) {
    # on log axes all x, xmin, xmax must be > 0
    dat <- dplyr::filter(dat, .x > 0, .xmin > 0, .xmax > 0)
  }
  
  # finite x-span for background stripes
  x_all <- c(dat$.x, dat$.xmin, dat$.xmax)
  x_all <- x_all[is.finite(x_all)]
  if (length(x_all) == 0) x_all <- c(1, 1)
  xr <- range(x_all, na.rm = TRUE)
  if (xr[1] == xr[2]) xr <- xr + c(-0.05, 0.05) * ifelse(xr[1] == 0, 1, abs(xr[1]))
  stripe_tmp <- transform(stripe_df, xmin = xr[1], xmax = xr[2])
  adjt<-1.3
  stripe_tmp$xmin <- xr[1] / adjt
  stripe_tmp$xmax <- xr[2] * adjt
  
  p <- ggplot(dat, aes(x = .data[[xvar]], y = Trait_ord)) +
    geom_rect(data = stripe_tmp,
              aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
              inherit.aes = FALSE, fill = "grey95") +
    geom_errorbarh(aes(xmin = .data[[xmin]], xmax = .data[[xmax]]),
                   height = 0.3, colour = point_col) +
    geom_point(size = 2.8, colour = point_col) +
    { if (length(ref_at) > 0) geom_vline(xintercept = ref_at, linetype = "dashed") else NULL } +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      axis.title.y     = element_blank(),
      axis.title.x     = element_blank(),
      axis.text.y      = if (show_y) element_text(size = 10) else element_blank(),
      axis.ticks.y     = if (show_y) element_line() else element_blank(),
      axis.text.x      = element_text(size = 10),
      axis.ticks.length= unit(4, "pt"),
      legend.position  = "none"
    )
  
  if (isTRUE(log_scale)) {
    # Data-aware log ticks using the same base as the annotation ticks
    p <- p +
      scale_x_continuous(
        trans        = scales::log_trans(base = log_base),
        breaks       = scales::log_breaks(base = log_base, n = tick_n),
        labels       = scales::label_number(),
        expand  = expansion(mult = 0) # keep plain numbers for ratios
        # (omit minor_breaks; default minors are fine with log_trans)
      ) +
      annotation_logticks(sides = "b", base = log_base,
                          short = unit(1.5, "pt"),     # length of minor ticks
                          mid   = unit(2.0, "pt"),     # length of mid ticks
                          long  = unit(3.0, "pt"))    # length of major ticks)  # ticks match the scale base
  } else {
    stripe_tmp <- transform(stripe_df, xmin = xr[1], xmax = xr[2])
    p <- ggplot(dat, aes(x = .data[[xvar]], y = Trait_ord)) +
      geom_rect(data = stripe_tmp,
                aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
                inherit.aes = FALSE, fill = "grey95") +
      geom_errorbarh(aes(xmin = .data[[xmin]], xmax = .data[[xmax]]),
                     height = 0.3, colour = point_col) +
      geom_point(size = 2.8, colour = point_col) +
      { if (length(ref_at) > 0) geom_vline(xintercept = ref_at, linetype = "dashed") else NULL } +
      coord_cartesian(clip = "off") +
      theme_classic() +
      theme(
        axis.title.y     = element_blank(),
        axis.title.x     = element_blank(),
        axis.text.y      = if (show_y) element_text(size = 10) else element_blank(),
        axis.ticks.y     = if (show_y) element_line() else element_blank(),
        axis.text.x      = element_text(size = 10),
        axis.ticks.length= unit(4, "pt"),
        legend.position  = "none"
      ) 
  }
  
  p
}

## generate each panel
p_h2 <- make_panel(h2_ME, "ratio_h2", "min_h2", "max_h2",
                   show_y = TRUE,  log_scale = TRUE, ref_at = 1)
p_me <- make_panel(h2_ME, "ratio_ME", "min_ME", "max_ME",
                   show_y = FALSE, log_scale = TRUE, log_base = 10, ref_at = 1) 
p_rg <- make_panel(h2_ME, "rg", "min_rg", "max_rg",
                   show_y = FALSE, log_scale = FALSE, ref_at = 1)  

##Combine
fig1 <- p_h2 + p_me + p_rg + plot_layout(widths = c(1, 1, 1))
#ggsave("fig1.png", fig1, width = 823, height = 211)

## Panel A with axis break (linear scale)
h2_lin <- h2_ME %>%
  mutate(
    se_h2_test =  ratio_h2 * sqrt( (h2_se_AoU / h2_AoU)^2 + (h2_se_UKB / h2_UKB)^2 ),
    lo = ratio_h2 - 1.96 * se_h2_test,
    hi = ratio_h2 + 1.96 * se_h2_test
  ) %>%
  arrange(desc(ratio_h2)) %>%
  mutate(Trait_f = factor(Trait, levels = Trait))

h2_lin  <- h2_lin  %>%
  filter(is.finite(ratio_h2), is.finite(lo), is.finite(hi)) %>%
  arrange(desc(ratio_h2)) %>%                    # larger ratios at bottom
  mutate(Trait_f = factor(Trait, levels = Trait))

t2d_name <- h2_lin$Trait[which.max(h2_lin$ratio_h2)]

# stripes
levs <- levels(h2_lin$Trait_f)
n <- length(levs)
stripe_df <- data.frame(
  y    = seq_len(n),
  ymin = seq_len(n) - 0.5,
  ymax = seq_len(n) + 0.5
) %>% filter(y %% 2 == 0)

pA_break<-ggplot(h2_lin, aes(x = ratio_h2, 
                             y = Trait_f)) +
  geom_rect(
    data = stripe_df,
    aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
    inherit.aes = FALSE,
    fill = "grey95"
  ) +
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 height = 0.30, linewidth = 0.5, colour = "#B56727") +
  geom_point(size = 2.6, colour = "#B56727") +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
  # axis break 
  scale_x_break(breaks = c(2.3, 4), scales = 0.25, space = 0.001) +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 10, 40),
                     expand = expansion(mult = c(0.02, 0.06))) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y      = element_blank(),
    axis.title.x      = element_blank(),
    axis.text.y       = element_text(size = 10),
    axis.ticks.length = unit(4, "pt"),
    legend.position   = "none",
    axis.title.x.top = element_blank(),
    axis.text.x.top  = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.x.top  = element_blank()
  )
pA_break

### Comparative deming reg
library(deming)

dr_h2<- deming(h2_AoU ~ 0 + h2_UKB, data = h2_ME,
              xstd = h2_ME$h2_se_UKB ,   # std. errors in X
              ystd = h2_ME$h2_se_AoU)   # std. errors in Y
summary(dr_h2)
coef(dr_h2)          # slope and intercept
1-coef(dr_h2)[[2]]



########## FIG 2

### load data 
aou_bins <-  fread("bin_all_aou.txt", sep="\t")
ukb_bins <-  fread("bin_all_ukb.txt", sep="\t")
fg_bins <-  fread("bin_all_fg.txt", sep="\t")

fg_threshold_001 <- fread("threshold_all_fg.txt")
fg_threshold_01 <- fread("threshold_all_fg_01.txt")
fg_threshold_1 <- fread("threshold_all_fg_1.txt")

ukb_threshold_001 <- fread("threshold_all_ukb.txt")
ukb_threshold_01 <- fread("threshold_all_ukb_01.txt")
ukb_threshold_1 <- fread("threshold_all_ukb_1.txt")

aou_threshold_001 <- fread("threshold_all_aou.txt")
aou_threshold_01 <- fread("threshold_all_aou_01.txt")
aou_threshold_1 <- fread("threshold_all_aou_1.txt")

### prep data
ukb_bins$trait <- gsub("standing_height", "height", ukb_bins$trait) #rename to height

ukb_binned <- ukb_bins %>% arrange(tolower(trait))
aou_binned <- aou_bins %>% arrange(tolower(trait))
fg_binned <- fg_bins %>% arrange(tolower(trait))

ukb_aou_traits <- c("Alzheimer's disease", "Asthma", "Basophil percentage", "BMI",
                    "Height","Monocyte percentage", "Neutrophil percentage", "Schizophrenia",
                    "Type 1 diabetes", "Type 2 diabetes", "Weight")

fg_traits <- c("Alzheimer's disease", "Asthma", 
               "Schizophrenia","Type 1 diabetes", 
               "Type 2 diabetes")

ukb_binned$trait<-rep(ukb_aou_traits, each = 22)
aou_binned$trait<-rep(ukb_aou_traits, each = 22)
fg_binned$trait<-rep(fg_traits, each = 22)

fg_threshold <- rbind(fg_threshold_001,fg_threshold_01,fg_threshold_1)
ukb_threshold <- rbind(ukb_threshold_001,ukb_threshold_01,ukb_threshold_1)
aou_threshold <- rbind(aou_threshold_001,aou_threshold_01,aou_threshold_1)

ukb_threshold$trait <- gsub("standing_height", "height", ukb_threshold$trait)

fg_threshold<- fg_threshold %>% arrange(tolower(trait))
ukb_threshold<- ukb_threshold %>% arrange(tolower(trait))
aou_threshold<- aou_threshold %>% arrange(tolower(trait))

fg_threshold$trait <- rep(fg_traits, each = 6)
ukb_threshold$trait <- rep(ukb_aou_traits, each = 6)
aou_threshold$trait <- rep(ukb_aou_traits, each = 6)

fg_threshold$database <- "FG"
ukb_threshold$database<- "UKB"
aou_threshold$database <- "AoU"

threshold_all <- rbind(ukb_threshold,aou_threshold,fg_threshold)


### Make plot functions

## set bins
bin_levels_maf <- c("bin1", "bin2", "bin3", "bin4", "bin5", "bin6", 
                    "bin7", "bin8", "bin9", "bin10", "bin11")

## function for Panel A
plot_sb <- function(ukb_dat, aou_dat, fg_dat=NULL,
                         trait_n, type_n, method_n,
                         bin_levels,
                         facet = FALSE,
                         bottom_pane = FALSE) {
  if (is_null(fg_dat)){
    comb_dat <- bind_rows(
      ukb_dat %>% mutate(databank = "UKB"),
      aou_dat %>% mutate(databank = "AoU")
    ) %>%
      filter(trait  == trait_n,
             type   == type_n,
             method == method_n) %>%
      mutate(bin = factor(bin, levels = bin_levels))
    
    gg <- ggplot(comb_dat,
                 aes(x = bin, y = mean_eta,
                     group = databank, colour = databank)) +
      geom_line() +
      geom_errorbar(aes(ymin = mean_eta - 1.96 * se,
                        ymax = mean_eta + 1.96 * se),
                    width = .3) +
      geom_point(size = 1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_colour_manual(values = c(UKB = "#D55E00", 
                                     AoU = "#009E73")) +
      ylim(-0.03,1) +
      scale_x_discrete(labels = comb_dat$bin_range) +
      theme_bw()
    
  } else {
    comb_dat <- bind_rows(
      ukb_dat %>% mutate(databank = "UKB"),
      aou_dat %>% mutate(databank = "AoU"),
      fg_dat %>% mutate(databank = "FG")
    ) %>%
      filter(trait  == trait_n,
             type   == type_n,
             method == method_n) %>%
      mutate(bin = factor(bin, levels = bin_levels))
    
    gg <- ggplot(comb_dat,
                 aes(x = bin, y = mean_eta,
                     group = databank, colour = databank)) +
      geom_line() +
      geom_errorbar(aes(ymin = mean_eta - 1.96 * se,
                        ymax = mean_eta + 1.96 * se),
                    width = .3) +
      geom_point(size = 1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_colour_manual(values = c(UKB = "#D55E00", 
                                     AoU = "#009E73", 
                                     FG = "#0072B2")) +
      ylim(-0.03,1) +
      scale_x_discrete(labels = comb_dat$bin_range) +
      theme_bw()
  }
  
  if (bottom_pane) {
    gg<- gg + theme(axis.text.x = element_text(size =  8, angle = 45, 
                                               hjust = 1),
                    axis.text.y = element_text(size = 8),
                    axis.title.y = element_blank(),
                    axis.title.x = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.title  = element_blank(),
                    panel.border = element_rect(colour = "black", 
                                                fill = NA, linewidth = 1),
                    legend.position = "none")
  } else {
    gg<- gg + theme(axis.text.x = element_blank(),
                    axis.text.y = element_text(size = 8),
                    axis.title.y = element_blank(),
                    axis.title.x = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.title  = element_blank(),
                    panel.border = element_rect(colour = "black", 
                                                fill = NA, linewidth = 1),
                    legend.position = "none")
  }
  
  
  if (facet) gg <- gg + facet_wrap(~databank, nrow = 1)
  
  gg
  
}

## plot panel A
scz_bin <- plot_sb(ukb_dat = ukb_binned,
                   aou_dat = aou_binned,
                   fg_dat = fg_binned,
                   trait_n = "Schizophrenia", 
                   type_n = "MAF", 
                   method_n = "sig", 
                   bin_levels = bin_levels_maf) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.5), "cm"))

baso_bin <- plot_sb(ukb_dat = ukb_binned,
                   aou_dat = aou_binned,
                   trait_n = "Basophil percentage", 
                   type_n = "MAF", 
                   method_n = "sig", 
                   bin_levels = bin_levels_maf) +
  theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm"))

t2d_bin <- plot_sb(ukb_dat = ukb_binned,
                   aou_dat = aou_binned,
                   fg_dat = fg_binned,
                   trait_n = "Type 2 diabetes", 
                   type_n = "MAF", 
                   method_n = "sig", 
                   bin_levels = bin_levels_maf,
                   bottom_pane = TRUE) +
  theme(plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm"))


# panel A grid
fig2a <- plot_grid(scz_bin, baso_bin , t2d_bin, 
          nrow = 3,
          rel_heights = c(1.2,1.1,1.75))

#ggsave("fig2_a.png", fig2a, width = 300, height = 435)


## Panel B

# plotting function
forest_t <- function(dt,tp,th,
                     higher=FALSE){
    p_dat <- dt %>% 
      filter(T_low == th & method==tp) %>%
      select(trait, mean_eta, lower_ci, upper_ci, database) %>%
      mutate(
        ## order traits by the *average* eta
        trait = reorder(trait, mean_eta, FUN = median)
    )
  
  dodge <- position_dodge(width = 0.6)          
  trait_levels <- levels(p_dat$trait)
  y_breaks <- seq(1.5, length(trait_levels) - 0.5, by = 1)
  
  if (higher==TRUE){
    sort_dat<-threshold_all[T_low == 0.001 & method==tp]
    p<-ggplot(p_dat,
              aes(x = mean_eta,
                  y = reorder(trait, sort_dat$mean_eta, FUN = median),
                  colour = database))
  } else {
    p<-ggplot(p_dat,
              aes(x = mean_eta,
                  y = trait,
                  colour = database))
  }

    p <- p + geom_point(size = 2.8, position = dodge) +
    geom_errorbarh(aes(xmin = lower_ci,
                       xmax = upper_ci),
                   height = 0.5,
                   position = dodge) +
    geom_hline(yintercept = y_breaks, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = c(UKB = "#D55E00", 
                                   AoU = "#009E73", 
                                   FG = "#0072B2")) +
    theme_classic() +
    theme(
      axis.title.y  = element_blank(),
      axis.title.x  = element_blank(),
      axis.text.y = element_text(size=10),
      axis.text.x = element_text(size=10),
      legend.title  = element_blank(),
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) 
  p
}

# plot panel B
fig2b<-forest_t(threshold_all,tp = "sig",th = 0.001)
fig2b

#ggsave("fig2_b.png", fig2b, width = 487, height = 521)

########## FIG 3 (simulation figure)
# see simulation script for plotting ()

########## FIG 4

### load data
moment_fg <- fread("fg_moments_1.csv", sep=",")
moment_aou <- fread("aou_moments_1.tsv", sep="\t")
moment_ukb <- fread("ukb_moments_1.csv", sep=",")

### prep data
moment_fg$trait <- gsub("F5_SCHZPHR", "schizophrenia", moment_fg$trait)
moment_fg$trait <- gsub("G6_ALZHEIMER", "alzheimers", moment_fg$trait)
moment_fg$trait <- gsub("J10_ASTHMA", "asthma", moment_fg$trait)

moment_fg <- moment_fg[c(-1,-2,-3),]

moment_fg  <- moment_fg  %>% arrange(tolower(trait))
moment_aou  <- moment_aou  %>% arrange(tolower(trait))
moment_ukb  <- moment_ukb  %>% arrange(tolower(trait))

moment_fg$trait  <- fg_traits
moment_aou$trait  <- ukb_aou_traits
moment_ukb$trait  <- ukb_aou_traits

fg_eta_skew <- merge(fg_threshold, moment_fg, by="trait")
ukb_eta_skew <- merge(ukb_threshold, moment_ukb, by="trait")
aou_eta_skew <- merge(aou_threshold, moment_aou, by="trait")

eta_skew <- bind_rows(fg_eta_skew, ukb_eta_skew, aou_eta_skew)

### Regression model

##  function
quad_model <- function(data){
  
  # fit pooled
  fit_pool <- lm(logit_p ~ poly(mu3_std, 2, raw = TRUE),
                 data = data)
  
  # add database intercepts
  fit_db   <- lm(logit_p ~ poly(mu3_std, 2, raw = TRUE) + database,
                 data = data)    
  
  # add database intercepts and slopes
  fit_db_int <- lm(logit_p ~ poly(mu3_std, 2, raw = TRUE) * database,
                   data = data)
  
  # summary of R2
  tbl_R2 <- tibble::tibble(
    model   = c("Pooled", "DB intercepts", "DB intercepts + slopes"),
    df      = c(fit_pool$df.residual,
                fit_db$df.residual,
                fit_db_int$df.residual),
    R2      = c(summary(fit_pool)$r.squared,
                summary(fit_db)$r.squared,
                summary(fit_db_int)$r.squared),
    adj_R2  = c(summary(fit_pool)$adj.r.squared,
                summary(fit_db)$adj.r.squared,
                summary(fit_db_int)$adj.r.squared)
  ) |>
    dplyr::mutate(
      delta_R2   = R2  - R2[1],        # gain over pooled model
      delta_adj  = adj_R2 - adj_R2[1]  # penalised gain
    )
  
  # correlation
  ct <- cor.test(data$mean_eta, 
                 data$mu3_std, 
                 method="spearman")
  
  #results object
  res<-list("fit_pooled" = fit_pool,
            "fit_db" = fit_db,
            "fit_db_int" = fit_db_int,
            "Fit pooled summary"=summary(fit_pool),
            "Fit db summary"=summary(fit_db),
            "Fit db+slopes summary"=summary(fit_db_int),
            "R2 table"=tbl_R2,
            "correlation"=ct)
  
  return(res)
}

### Plotting function
plot_m <- function(dat,mod, ribbon=TRUE){
  # make grid
  grid <- data.frame(
    mu3_std = seq(min(dat$mu3_std), max(dat$mu3_std), length.out = 400)
  )
  
  # get predicted values and SEs on logit scale
  pred <- predict(mod$fit_pooled, newdata = grid, se.fit = TRUE)
  grid <- grid %>%
    mutate(
      z_hat   = pred$fit,
      z_low   = z_hat - 1.96 * pred$se.fit,
      z_high  = z_hat + 1.96 * pred$se.fit,
      
      p_hat   = plogis(z_hat),
      p_low   = plogis(z_low),
      p_high  = plogis(z_high),
      
      eta_hat = 2 * p_hat  - 1,
      eta_low = 2 * p_low  - 1,
      eta_high= 2 * p_high - 1
    )
  
  #plot
  
  p<-ggplot(dat,
         aes(x = mu3_std, y = mean_eta,
             colour = database, label = trait)) +
    geom_point(size = 2) +
    geom_text_repel(size = 3, max.overlaps = 20, 
                    show.legend = FALSE) +
    scale_colour_manual(values = c(UKB = "#D55E00", 
                                   AoU = "#009E73", 
                                   FG = "#0072B2")) +
    geom_line(
      data = grid,
      aes(mu3_std, eta_hat),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.5
    ) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title  = element_blank(),
          legend.position = "none")
  
  if(ribbon==TRUE) {
    p<- p + geom_ribbon(
      data = grid,
      aes(x = mu3_std, ymin = eta_low, ymax = eta_high),
      inherit.aes = FALSE,
      fill = "grey80", alpha = 0.4
    ) 
  }
  
  return(p)
}

## convert to probabilities 
logit_func <- function(p, eps = 1e-6) {
  qlogis(pmin(pmax(p, eps), 1 - eps))
}

eta_skew <- eta_skew %>% 
  mutate(p = (mean_eta + 1) / 2) %>% 
  mutate(logit_p = logit_func(p))  

eta_skew$database <- factor(eta_skew$database,
                                levels = c("UKB", "AoU", "FG"))

## fit pooled (lowest p-value SNPs, MAF threshold T_low <= 0.001)
eta_skew_sig_001 <- eta_skew[method=="sig" & T_low == 0.001]
model_sig_001 <- quad_model(eta_skew_sig_001)

print(model_sig_001$`Fit pooled summary`)
print(model_sig_001$`Fit db summary`)
print(model_sig_001$`Fit db+slopes summary`)
print(model_sig_001$`R2 table`)
print(model_sig_001$correlation)

## fit pooled (random SNPs, MAF threshold T_low <= 0.001)
eta_skew_rand_001 <- eta_skew[method=="random" & T_low == 0.001]
model_rand_001 <- quad_model(eta_skew_rand_001)

print(model_rand_001$`Fit pooled summary`)
print(model_rand_001$`Fit db summary`)
print(model_rand_001$`Fit db+slopes summary`)
print(model_rand_001$`R2 table`)
print(model_rand_001$correlation)

# plot
sig_eta <- plot_m(eta_skew_sig_001, model_sig_001, 
                          ribbon=FALSE)
rand_eta <- plot_m(eta_skew_rand_001, 
                   model_rand_001,
                   ribbon=FALSE) + 
  ylim(-0.1,1) + 
  xlim(0,60)

fig4 <- plot_grid(sig_eta, rand_eta)

#ggsave("fig4.png", fig4, width = 570, height = 376)

##################################################
## Supplemental figures


########## FIG S1 (UKB and AoU quantitative distributions)

### load data
aou_bmi <- fread("aou_bmi.tsv", sep="\t")
aou_height <- fread("aou_height.tsv", sep="\t")
aou_weight <- fread("aou_weight.tsv", sep="\t")
aou_baso <- fread("aou_basophil_percentage.tsv", sep="\t")
aou_mono <- fread("aou_monocyte_percentage.tsv", sep="\t")
aou_neut <- fread("aou_neutrophil_percentage.tsv", sep="\t")

ukb_phenos <- fread("ukb_phenotypes.tsv", sep="\t")

### prep data
aou_list <- list(aou_height, aou_weight, aou_bmi,
                 aou_baso, aou_neut, aou_mono)
aou_list<- lapply(aou_list, function(x) { x$person_id <- as.integer(x$person_id); x })
aou_phenos <- Reduce(function(x, y) merge(x, y, by = "person_id", all = TRUE), aou_list)
names(aou_phenos) <- c("eid","height","weight","bmi","baso","neut","mono")

cols <- c("eid", "Height (cm)","Weight (kg)","BMI","Basophil percentage (median)",
          "Neutrophil percentage (median)","Monocyte percentage (median)")

names(aou_phenos) <- cols
names(ukb_phenos) <- cols

### Plot S1 
trait_names <- cols[-1]

aou_long <- aou_phenos %>%
  pivot_longer(-eid, names_to = "trait", values_to = "value") %>%
  mutate(source = "AoU")

ukb_long <- ukb_phenos %>%
  pivot_longer(-eid, names_to = "trait", values_to = "value") %>%
  mutate(source = "UKB") 

pheno_plot <- bind_rows(aou_long, ukb_long) %>%
  mutate(
    trait  = factor(trait, levels = trait_names),
    source = factor(source, levels = c("UKB","AoU")),
    value  = as.numeric(value)
  )

fig_s1 <- ggplot(pheno_plot, aes(x = value, 
                            fill = source, color = source)) +
  geom_density(alpha = 0.30, size = 0.6, adjust = 1.5, 
               position = "identity", na.rm = TRUE) +
  facet_wrap(~ trait, ncol = 3, scales = "free") +
  scale_fill_manual(values = c(AoU = "#007D58", UKB = "#B91C1C")) +
  scale_color_manual(values = c(AoU = "#009E73", UKB = "#D55E00")) +
  labs(x = NULL, y = NULL,
       title = NULL,
       fill = NULL, color = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank()
  )

fig_s1
# ggsave("aou_ukb_distributions.png", fig_s1, width = 1246, height = 741)



########## FIG S2 (Sign bias vs MAF bins, lowest p SNPs)

### plotting function
plot_sb_all <- function(ukb_dat, aou_dat, fg_dat,
                          type_n, method_n,
                          bin_levels) {
  comb_dat <- bind_rows (
    ukb_dat %>% mutate(Database = "UKB"),
    aou_dat %>% mutate(Database = "AoU"),
    fg_dat %>% mutate(Database = "FinnGen")
  ) %>%
    filter(type   == type_n,
           method == method_n) %>%
    mutate(bin = factor(bin, levels = bin_levels))
  
  gg <- ggplot(comb_dat,
               aes(x = bin, y = mean_eta,
                   group = Database, colour = Database)) +
    geom_line() +
    geom_errorbar(aes(ymin = mean_eta - 1.96 * se,
                      ymax = mean_eta + 1.96 * se),
                  width = .3) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = c(UKB = "#D55E00", 
                                   AoU = "#009E73", 
                                   FinnGen = "#0072B2")) +
    ylim(-0.18,1) +
    ylab("Sign bias") + 
    xlab("Minor allele frequency bin") +
    scale_x_discrete(labels = comb_dat$bin_range) +
    theme_classic() 
  
  
  gg <- gg + facet_wrap(~trait, nrow = 4) +
    theme(panel.grid.major = element_line(colour = "grey90"),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
          axis.text.x = element_text(size = 8, angle = 45, 
                                     hjust = 1),
          panel.spacing.x = unit(5, "mm"))
  
  gg
}

### Plot figure
fig_s2 <- plot_sb_all(ukb_dat = ukb_binned,
                      aou_dat = aou_binned,
                      fg_dat = fg_binned,
                      type_n = "MAF", 
                      method_n = "sig", 
                      bin_levels = bin_levels_maf)

# ggsave("fig_s2.png", fig_s2, width = 956, height = 640)




########## FIG S3 (Sign bias vs MAF bins, randomly sampled SNPs)

### Plot figure
fig_s3 <- plot_sb_all(ukb_dat = ukb_binned,
                      aou_dat = aou_binned,
                      fg_dat = fg_binned,
                      type_n = "MAF", 
                      method_n = "random", 
                      bin_levels = bin_levels_maf)

# ggsave("fig_s3.png", fig_s3, width = 956, height = 640)



########## FIG S4 (Sign bias at thresholds T_low 0.001, 0.01, 0.1; lowest p SNPs)

# plots, 3 thresholds

T_001 <- forest_t(threshold_all, 
         tp = "sig",
         th = 0.001)

T_01 <- forest_t(threshold_all, 
         tp = "sig",
         th = 0.01,
         higher = TRUE) +
  theme(axis.text.y = element_blank())

T_1 <- forest_t(threshold_all, 
         tp = "sig",
         th = 0.1,
         higher = TRUE) +
  theme(axis.text.y = element_blank())

# assemble plot grid
fig_s4 <- plot_grid(T_001, T_01, T_1, 
                   ncol = 3,
                   rel_widths = c(1,0.7,0.7))

# ggsave("fig_s4.png", fig_s4, width = 958, height = 332)



########## FIG S5 (Sign bias at thresholds T_low 0.001, 0.01, 0.1; randomly sampled SNPs)

# plots, 3 thresholds

T_001_r <- forest_t(threshold_all, 
                  tp = "random",
                  th = 0.001)

T_01_r <- forest_t(threshold_all, 
                 tp = "random",
                 th = 0.01,
                 higher = TRUE) +
  theme(axis.text.y = element_blank())

T_1_r <- forest_t(threshold_all, 
                tp = "random",
                th = 0.1,
                higher = TRUE) +
  theme(axis.text.y = element_blank())

# assemble plot grid
fig_s5 <- plot_grid(T_001_r, T_01_r, T_1_r, 
                    ncol = 3,
                    rel_widths = c(1,0.7,0.7))

# ggsave("fig_s5.png", fig_s5, width = 958, height = 332)



########## FIG S6 (Sign bias vs skew in UKB at thresholds T_low 0.001, 0.01, 0.1, lowest p)

### prep data
eta_skew_sig_01 <- eta_skew[method=="sig" & T_low == 0.01]
eta_skew_sig_1 <- eta_skew[method=="sig" & T_low == 0.1]

### fit pooled (random sampling, three thresholds)
model_sig_01 <- quad_model(eta_skew_sig_01)
model_sig_1 <- quad_model(eta_skew_sig_1)

### plot
sig_001 <- plot_m(eta_skew_sig_001, model_sig_001) + ylim(-0.1,1) + xlim(0,60)
sig_01 <- plot_m(eta_skew_sig_01, model_sig_01)+ ylim(-0.1,1) + xlim(0,60)
sig_1 <- plot_m(eta_skew_sig_1, model_sig_1)+ ylim(-0.1,1)+ xlim(0,60)

# assemble plot grid
fig_s6 <- plot_grid(sig_001, sig_01, sig_1, 
                    ncol = 3)

# ggsave("fig_s6.png", fig_s6, width = 1630, height = 421)



########## FIG S7 (Sign bias vs skew in UKB at thresholds T_low 0.001, 0.01, 0.1, random sampling)

### prep data
eta_skew_rand_001 <- eta_skew[method=="random" & T_low == 0.001]
eta_skew_rand_01 <- eta_skew[method=="random" & T_low == 0.01]
eta_skew_rand_1 <- eta_skew[method=="random" & T_low == 0.1]

### fit pooled (random sampling, three thresholds)
model_rand_001 <- quad_model(eta_skew_rand_001)
model_rand_01 <- quad_model(eta_skew_rand_01)
model_rand_1 <- quad_model(eta_skew_rand_1)

### plot
rand_001 <- plot_m(eta_skew_rand_001, model_rand_001) + ylim(-0.1,1) + xlim(0,60)
rand_01 <- plot_m(eta_skew_rand_01, model_rand_01)+ ylim(-0.1,1) + xlim(0,60)
rand_1 <- plot_m(eta_skew_rand_1, model_rand_1)+ ylim(-0.1,1)+ xlim(0,60)

# assemble plot grid
fig_s7 <- plot_grid(rand_001, rand_01, rand_1, 
                    ncol = 3)

# ggsave("fig_s7.png", fig_s7, width = 1630, height = 421)



########## FIG S8 (Results comparing IRNT traits to raw traits)

#load data
ukb_threshold_raw <- fread("threshold_all_ukb.txt")
ukb_threshold_irnt <- fread("irnt/threshold_ukb_irnt.txt")

ukb_threshold_raw$trait <- gsub("neutrophill_percentage", "neutrophil_percentage", ukb_threshold_raw$trait)

ukb_threshold_raw  <- ukb_threshold_raw[method == "sig"]
ukb_threshold_raw  <- ukb_threshold_raw[c(3,4,5,6,8,11),]

ukb_threshold_raw$datatype <- "Raw"
ukb_threshold_irnt$datatype <- "IRN"

#bind
threshold_both <- bind_rows(ukb_threshold_raw, ukb_threshold_irnt)

#plot threshold
threshold_both$trait<- c("Basophil percentage", "BMI",  "Monocyte percentage",  
                         "Neutrophil percentage", "Height", "Weight",
                         "Basophil percentage", "BMI", "Monocyte percentage",
                         "Neutrophil percentage", "Height", "Weight")  

plot_dat <- threshold_both %>%            
  select(trait, mean_eta, lower_ci, upper_ci, datatype) %>%
  mutate(
    trait = reorder(trait, mean_eta, FUN = median)
  )

trait_levels <- levels(as.factor(plot_dat$trait))
y_breaks <- seq(1.5, length(trait_levels) - 0.5, by = 1)

#w: 487, h:521
ggplot(plot_dat,
       aes(x = mean_eta,
           y = trait,
           colour = datatype)) +
  geom_point(size = 2.8) +
  geom_errorbarh(aes(xmin = lower_ci,
                     xmax = upper_ci),
                 height = 0.2) +
  geom_hline(yintercept = y_breaks, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = c(Raw = "#000000", IRN = "darkred")) +
  theme_classic() +
  theme(
    axis.title.y  = element_blank(),
    axis.title.x  = element_blank(),
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(size=10),
    legend.title  = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )






#load maf bin data
ukb_bin_raw <- fread("../bin_ukb_quant.txt")
ukb_bin_irnt <- fread("irnt/bin_ukb_irnt.txt")

ukb_bin_raw$trait <- gsub("neutrophill_percentage", "neutrophil_percentage", ukb_bin_raw$trait)

bin_levels_maf <- c("bin1", "bin2", "bin3", "bin4", "bin5", "bin6", 
                    "bin7", "bin8", "bin9", "bin10", "bin11")

#bin plotting function
plot_sb_dual_irnt <- function(ukb_dat, irnt_dat,
                              trait_a, trait_n, type_n, method_n,
                              bin_levels,
                              cols = c(Raw = "#000000", IRN = "darkred"),
                              facet = FALSE,
                              bottom_pane = FALSE) {
  
  comb_dat <- bind_rows(
    ukb_dat %>% mutate(datatype = "Raw"),
    irnt_dat %>% mutate(datatype = "IRN"),
  ) %>%
    filter(trait  == trait_a,
           type   == type_n,
           method == method_n) %>%
    mutate(bin = factor(bin, levels = bin_levels))
  
  gg <- ggplot(comb_dat,
               aes(x = bin, y = mean_eta,
                   group = datatype, colour = datatype)) +
    geom_line() +
    geom_errorbar(aes(ymin = mean_eta - 1.96 * se,
                      ymax = mean_eta + 1.96 * se),
                  width = .3) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = cols) +
    scale_x_discrete(labels = comb_dat$bin_range) +
    ggtitle(trait_n) +
    ylim(-0.15,1) +
    xlab("MAF bin") + ylab("Sign bia") +
    theme_bw()
  
  gg<- gg + theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
                  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                  legend.position  = "none")
  
  
  if (facet) gg <- gg + facet_wrap(~databank, nrow = 1)
  
  gg
}

#### bin plots
bmi <- plot_sb_dual_irnt(ukb_dat = ukb_bin_raw,
                         irnt_dat = ukb_bin_irnt,
                         trait_a = "BMI",
                         trait_n = "BMI", 
                         type_n = "MAF",
                         method_n = "sig",
                         bin_levels = bin_levels_maf,
                         bottom_pane = FALSE) 

height <- plot_sb_dual_irnt(ukb_dat = ukb_bin_raw,
                            irnt_dat = ukb_bin_irnt,
                            trait_a = "standing_height",
                            trait_n = "Standing height", 
                            type_n = "MAF",
                            method_n = "sig",
                            bin_levels = bin_levels_maf,
                            bottom_pane = FALSE) 

weight <- plot_sb_dual_irnt(ukb_dat = ukb_bin_raw,
                            irnt_dat = ukb_bin_irnt,
                            trait_a = "weight",
                            trait_n = "Weight", 
                            type_n = "MAF",
                            method_n = "sig",
                            bin_levels = bin_levels_maf,
                            bottom_pane = FALSE) 

mono <- plot_sb_dual_irnt(ukb_dat = ukb_bin_raw,
                          irnt_dat = ukb_bin_irnt,
                          trait_a = "monocyte_percentage",
                          trait_n = "Monocyte percentage", 
                          type_n = "MAF",
                          method_n = "sig",
                          bin_levels = bin_levels_maf,
                          bottom_pane = FALSE) 

baso <- plot_sb_dual_irnt(ukb_dat = ukb_bin_raw,
                          irnt_dat = ukb_bin_irnt,
                          trait_a = "basophil_percentage",
                          trait_n = "Basophil percentage", 
                          type_n = "MAF",
                          method_n = "sig",
                          bin_levels = bin_levels_maf,
                          bottom_pane = FALSE) 

neut <- plot_sb_dual_irnt(ukb_dat = ukb_bin_raw,
                          irnt_dat = ukb_bin_irnt,
                          trait_a = "neutrophil_percentage",
                          trait_n = "Neutrophil percentage", 
                          type_n = "MAF",
                          method_n = "sig",
                          bin_levels = bin_levels_maf,
                          bottom_pane = FALSE) 

#panels
plot_grid(height, bmi, weight,
          mono, baso, neut,
          ncol = 3, nrow = 2)
#w: 703, h:510



########## FIG S9 (Plotting genetic correlations for subsamples of UKB)

# read rg results
rg <- fread("../../null_rg_analysis.csv", sep=",") 
rg[, trait := factor(trait,
                     levels = c("Height", "Weight", "BMI",
                                "Monocyte percentage", "Basophil percentage", "Neutrophil percentage",
                                "Asthma", "Type 1 diabetes", "Type 2 diabetes", "Schizophrenia"))]
rg[, comp := factor(
  comp,
  levels = c("AOU_Neale", "UKB1_vs_UKB2", "UKB1_vs_AOU", "UKB2_vs_AOU"),
  labels = c("AoU vs Neale",
             "UKB1 vs UKB2",
             "UKB1 vs AoU",
             "UKB2 vs AoU")
)]

quant_traits  <- c("Height", "Weight", "BMI",
                   "Monocyte percentage", "Basophil percentage", "Neutrophil percentage")
binary_traits <- c("Asthma", "Type 1 diabetes", "Type 2 diabetes", "Schizophrenia")

rg_quant <- rg[trait %in% quant_traits]
rg_bin   <- rg[trait %in% binary_traits]

##Plot quant
p_quant <- ggplot(rg_quant, aes(x = comp, y = rg, color = comp)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, position = position_dodge(width = 0.4)) +
  geom_errorbar(
    aes(ymin = rg - 1.96 * SE,
        ymax = rg + 1.96 * SE),
    width = 0.2,
    position = position_dodge(width = 0.4)
  ) +
  facet_wrap(~ trait, ncol = 3) +
  coord_cartesian(ylim = c(0.01, 1.1)) +
  scale_color_discrete(name = "Comparison") +
  labs(
    x = "",
    y = "Genetic correlation (rg)",
    title = "Quantitative traits"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

p_quant

##plot binary
p_bin <- ggplot(rg_bin, aes(x = comp, y = rg, color = comp)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, position = position_dodge(width = 0.4)) +
  geom_errorbar(
    aes(ymin = rg - 1.96 * SE,
        ymax = rg + 1.96 * SE),
    width = 0.2,
    position = position_dodge(width = 0.4)
  ) +
  facet_wrap(~ trait, ncol = 2) +
  coord_cartesian(ylim = c(-1.0, 2.5)) +  # wider range for noisy binary traits
  scale_color_discrete(name = "Comparison") +
  labs(
    x = "",
    y = "Genetic correlation (rg)",
    title = "Binary traits"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

p_bin

plot_grid(p_quant,p_bin, nrow=1,
          rel_widths = c(1,0.7))

#w: 835, h: 417
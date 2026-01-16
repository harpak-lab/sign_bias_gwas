############################################################
## Plot simulation results
##
## Inputs:
##   1) sample RDS produced by --single:
##   2) out_reps TSV produced by --sim (plain text)
##
## Produces plots
##   Figure 3, pA, pB, pC, pD, pE
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
})

sample_rds_path <- "sim1_sample.rds"
out_reps_path   <- "sim1_out_reps.txt"  

# Load inputs
sims_run_f3 <- readRDS(sample_rds_path)
if (!is.list(sims_run_f3) || length(sims_run_f3) < 5) {
  stop("sample RDS doesn't look like sims_run_f3 <- list(Y_pop, Y_demo, G_inc, G_dec, gwas_cohort)")
}

Y_pop       <- sims_run_f3[[1]]
Y_demo      <- sims_run_f3[[2]]
G_inc       <- sims_run_f3[[3]]
G_dec       <- sims_run_f3[[4]]
gwas_cohort <- sims_run_f3[[5]]
beta <- sims_run_f3$pb %||% sims_run_f3[[6]]

out_reps <- read.delim(out_reps_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
out_reps <- as_tibble(out_reps)


# Panel A (pop + cohort density)
dens_pop <- density(Y_pop, adjust = 3)
dd_pop <- data.frame(x = dens_pop$x, y = dens_pop$y)

pA1 <- ggplot(dd_pop, aes(x, y)) +
  geom_area(fill = "#FDE992", color = NA) +
  geom_line(color = "black", linewidth = 0.6) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = 12)
  )

dens_samp <- density(Y_demo, adjust = 3, from = -5, to = 5)
dd_samp <- data.frame(x = dens_samp$x, y = dens_samp$y)
q90_samp <- as.numeric(quantile(Y_demo, 0.90))

pA2 <- ggplot(dd_samp, aes(x, y)) +
  geom_area(fill = "#EFEFEF", color = NA) +
  geom_area(data = subset(dd_samp, x >= q90_samp), fill = "orange3", color = NA) +
  geom_line(color = "black", linewidth = 0.6) +
  geom_vline(xintercept = q90_samp, linetype = "dashed") +
  xlim(-5, 5) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )

pA <- plot_grid(pA1, pA2, nrow = 2, align = "v", rel_heights = c(0.9, 1))


# Panel B (scatter Y ~ G_inc)
q90 <- as.numeric(stats::quantile(Y_demo, 0.90))
tail_flag <- Y_demo > q90

mk_panel2 <- function(G) {
  df <- data.frame(
    G = G,
    Y = Y_demo,
    tail = factor(ifelse(tail_flag, "Right tail", "Bulk"),
                  levels = c("Bulk", "Right tail"))
  )
  
  df$lev <- stats::hatvalues(stats::lm(Y ~ G, data = df))
  
  set.seed(1)
  df <- df[sample.int(nrow(df), min(nrow(df), 2500)), ]
  
  ggplot(df, aes(G, Y)) +
    geom_point(aes(size = lev, fill = tail),
               shape = 21, alpha = 0.8, stroke = 0.25, color = "gray20") +
    scale_fill_manual(values = c("Bulk" = "#BDBDBD", "Right tail" = "orange3"),
                      guide = "none") +
    scale_size_continuous(range = c(1.6, 4.6), guide = "none") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
    scale_x_continuous(breaks = 0:2, limits = c(-0.15, 2.15)) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.position = "none")
}

pB <- mk_panel2(G_inc)


# Panel C (carrier rates, bulk vs tail)
car_df <- tibble(
  snp   = rep(c("Trait decreasing", "Trait increasing"), each = 2),
  group = rep(c("Bulk", "Right tail"), times = 2),
  carrier_rate = c(
    mean((G_dec > 0)[!tail_flag]),
    mean((G_dec > 0)[ tail_flag]),
    mean((G_inc > 0)[!tail_flag]),
    mean((G_inc > 0)[ tail_flag])
  )
) |>
  mutate(
    snp   = factor(snp,   levels = c("Trait decreasing", "Trait increasing")),
    group = factor(group, levels = c("Bulk", "Right tail"))
  )

pC <- ggplot(car_df, aes(x = snp, y = carrier_rate, fill = group)) +
  geom_col(position = position_dodge(width = .6),
           width = .55, color = "gray20") +
  scale_fill_manual(values = c("Bulk" = "#BDBDBD",
                               "Right tail" = "orange3")) +
  geom_text(aes(label = sprintf("%.2f", carrier_rate)),
            position = position_dodge(width = .6),
            vjust = -0.35, size = 4) +
  geom_hline(yintercept = 0.51,
             linetype = "dashed",
             color = "gray40",
             linewidth = 0.6) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")


# Panel D 
if (!("sign_true" %in% names(gwas_cohort))) {
  s <- sign(beta[gwas_cohort$snp])
  s_minor <- ifelse(gwas_cohort$af <= 0.5, s, -s)
  gwas_cohort$sign_true <- factor(
    ifelse(s_minor > 0, "Trait increasing", "Trait decreasing"),
    levels = c("Trait decreasing", "Trait increasing")
  )
}
gwas_forD <- gwas_cohort |>
  mutate(sign_true = factor(sign_true, levels = c("Trait decreasing", "Trait increasing"))) |>
  filter(is.finite(se), is.finite(sign_true))

pD <- ggplot(gwas_forD, aes(x = sign_true, y = se)) +
  stat_summary(fun = mean, linewidth = 0.6) +
  stat_summary(fun.data = mean_cl_normal,
               geom = "errorbar", width = .3, linewidth = 0.6) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")


# Panel E (binned skew vs estimated/true sign bias)
need_cols <- c("skew_obs", "bias_sig", "true_bias_sample")
miss <- setdiff(need_cols, names(out_reps))
if (length(miss)) stop("out_reps missing required columns: ", paste(miss, collapse = ", "))

rng  <- range(out_reps$skew_obs, na.rm = TRUE)
brks <- seq(rng[1], rng[2], length.out = 9)

sum4 <- out_reps %>%
  mutate(bin = cut(skew_obs, breaks = brks, include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(
    x = mean(skew_obs, na.rm = TRUE),
    n = dplyr::n(),
    mean_est = mean(bias_sig, na.rm = TRUE),
    se_est   = sd(bias_sig,   na.rm = TRUE) / sqrt(pmax(n, 1)),
    mean_true = mean(true_bias_sample, na.rm = TRUE),
    se_true   = sd(true_bias_sample,   na.rm = TRUE) / sqrt(pmax(n, 1)),
    .groups = "drop"
  )

pE_plot <- bind_rows(
  sum4 %>% transmute(x, series = "Cohort true sign bias",
                     mean = mean_true, se = se_true, shape = "true"),
  sum4 %>% transmute(x, series = "Estimated sign bias",
                     mean = mean_est,  se = se_est,  shape = "est")
) %>%
  mutate(ymin = mean - 2*se, ymax = mean + 2*se)

pE <- ggplot(pE_plot, aes(x = x, y = mean, color = series, group = series)) +
  geom_hline(yintercept = 0, linewidth = 1, color = "#D5B60A") +
  geom_line(linewidth = 0.7) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                width = 0.035, linewidth = 0.5) +
  geom_point(aes(shape = shape), size = 3.0, stroke = 0.9) +
  scale_shape_manual(values = c(true = 17, est = 16), guide = "none") +
  scale_color_manual(values = c(
    "Cohort true sign bias" = "purple4",
    "Estimated sign bias"   = "black"
  )) +
  scale_x_continuous(breaks = 0:4, limits = c(0, 4)) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_text(size = 12),
    legend.position = "none"
  )

#Display plots
pA
pB
pC
pD
pE

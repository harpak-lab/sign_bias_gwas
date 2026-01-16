#!/usr/bin/env Rscript
##### Run ASH on GWAS sumstats

library("dplyr")
library("dtplyr")
library("tidyr")
library("data.table")
library("ashr")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_ash.R <cleaned_gwas_file> <out_file>\nExample: Rscript run_ash.R cleaned_gwas/asthma.cleaned.txt ash_out/asthma.minor_results.txt")
}

in_file  <- args[1]
out_file <- args[2]

cleaned_dat <- fread(in_file)

ash_minor_out <- ash(cleaned_dat$m_Beta, cleaned_dat$SE)
ash_results_minor <- ash_minor_out$result

minor_results <- cbind(cleaned_dat, ash_results_minor)
minor_results$eta <- minor_results$PositiveProb - minor_results$NegativeProb

fwrite(
  minor_results,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  showProgress = TRUE
)

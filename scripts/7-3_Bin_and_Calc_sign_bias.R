#!/usr/bin/env Rscript

######## Binning/threshold 

###############################################################################
# A single script implementing:
#  (1) Bin-based sign bias (MAF=11 bins)
#  (2) Threshold-based sign bias (low or/and high cutoff),
# for either "random" SNP selection or "sig" (lowest p-value) SNP selection.
#

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (!("--file" %in% args)) {
  stop("Please specify --file <path>.")
}
if (!("--type" %in% args)) {
  stop("Please specify --type MAF or DAF.")
}
if (!("--method" %in% args)) {
  stop("Please specify --method random or sig.")
}
if (!("--mode" %in% args)) {
  stop("Please specify --mode bin or threshold.")
}

file_arg   <- args[which(args=="--file")+1]
type_arg   <- args[which(args=="--type")+1]    # "MAF" or "DAF"
method_arg <- args[which(args=="--method")+1]  # "random" or "sig"
mode_arg   <- args[which(args=="--mode")+1]    # "bin" or "threshold"

# --low and --high for threshold approach
low_idx  <- which(args=="--low")
high_idx <- which(args=="--high")

T_low  <- if (length(low_idx)>0)  as.numeric(args[low_idx+1]) else NA
T_high <- if (length(high_idx)>0) as.numeric(args[high_idx+1]) else NA

# Trait name from file
trait_name <- sub(".*\\.(.*?)\\..*", "\\1", basename(file_arg))

cat("Running sign-bias script on:\n",
    "  File:   ", file_arg, "\n",
    "  Type:   ", type_arg, "\n",
    "  Method: ", method_arg, "\n",
    "  Mode:   ", mode_arg, "\n",
    "  Trait:  ", trait_name, "\n",
    "  T_low=  ", T_low, ", T_high= ", T_high, "\n")

suppressMessages(library(data.table))
suppressMessages(library(dplyr))

###############################################################################
# 1) Define bin breaks (UKB)
###############################################################################
bin_breaks_maf <- c(
  0.0, 5e-4, 1e-3,
  0.002, 0.003, 0.005, 0.01,
  0.02, 0.05, 0.1, 0.3, 0.5
)
bin_breaks_daf <- c(
  0.0, 5e-4, 1e-3,
  0.002, 0.003, 0.005, 0.01, 0.02, 0.05,
  0.1, 0.3, 0.5, 0.7, 0.9,
  0.95, 0.98, 0.99, 0.995, 0.998, 1.0
)

# bin_breaks_maf <- c(
#   0.0, 0.001, 0.002, 0.003, 
#   0.004, 0.005, 0.01, 0.02, 0.05, 
#   0.1, 0.3, 0.5                
# )
# 
# bin_breaks_daf <- c(
#   0.0, 0.001, 0.002, 0.003, 0.004, 0.005,
#   0.01, 0.02, 0.05,
#   0.1, 0.3, 0.5, 0.7, 0.9,
#   0.95, 0.98, 0.99, 0.995, 0.998, 1.0
# )

###############################################################################
# 2) Functions for BIN mode
###############################################################################

calc_bin_stats_random <- function(df, bin_breaks, N=1000, min_snps=10) {
  bin_labels <- paste0("bin", seq_len(length(bin_breaks)-1))
  df[, bin_id := cut(af, breaks=bin_breaks, labels=bin_labels, include.lowest=TRUE)]
  
  B <- length(bin_labels)
  bin_results <- vector("list", B)
  
  for (b_idx in seq_len(B)) {
    this_bin <- bin_labels[b_idx]
    sub_df   <- df[bin_id == this_bin]
    
    if (nrow(sub_df)==0) {
      bin_results[[b_idx]] <- data.frame(
        bin       = this_bin,
        bin_range = paste0("[", bin_breaks[b_idx], ",", bin_breaks[b_idx+1], ")"),
        n_snps    = 0,
        n_blocks  = 0,
        mean_eta  = NA_real_,
        se        = NA_real_,
        lower_ci  = NA_real_,
        upper_ci  = NA_real_
      )
      next
    }
    # Filter blocks with >= min_snps
    counts_block <- sub_df[,.N, by=block_id]
    valid_blocks <- counts_block[N>=min_snps]$block_id
    sub_df <- sub_df[block_id %in% valid_blocks]
    
    if (nrow(sub_df)==0) {
      bin_results[[b_idx]] <- data.frame(
        bin       = this_bin,
        bin_range = paste0("[", bin_breaks[b_idx], ",", bin_breaks[b_idx+1], ")"),
        n_snps    = 0,
        n_blocks  = 0,
        mean_eta  = NA_real_,
        se        = NA_real_,
        lower_ci  = NA_real_,
        upper_ci  = NA_real_
      )
      next
    }
    
    blocks_in_bin <- unique(sub_df$block_id)
    n_snps_bin    <- nrow(sub_df)
    n_blocks_bin  <- length(blocks_in_bin)
    
    block_map <- split(sub_df$eta, sub_df$block_id)
    
    replicate_etas <- numeric(N)
    for (i in seq_len(N)) {
      chosen <- sapply(block_map, function(etas) sample(etas,1))
      replicate_etas[i] <- sum(chosen)/sum(abs(chosen))
    }
    mean_eta <- mean(replicate_etas)
    sd_eta   <- sd(replicate_etas)
    l_ci     <- mean_eta - 1.96*sd_eta
    u_ci     <- mean_eta + 1.96*sd_eta
    
    bin_results[[b_idx]] <- data.frame(
      bin       = this_bin,
      bin_range = paste0("[",bin_breaks[b_idx],",",bin_breaks[b_idx+1],")"),
      n_snps    = n_snps_bin,
      n_blocks  = n_blocks_bin,
      mean_eta  = mean_eta,
      se        = sd_eta,
      lower_ci  = l_ci,
      upper_ci  = u_ci
    )
  }
  
  bin_df <- do.call(rbind, bin_results)
  rownames(bin_df) <- NULL
  return(bin_df)
}


calc_bin_stats_sig <- function(df, bin_breaks, n_boot=1000, min_snps=10) {
  bin_labels <- paste0("bin", seq_len(length(bin_breaks)-1))
  df[, bin_id := cut(af, breaks=bin_breaks, labels=bin_labels, include.lowest=TRUE)]
  
  B <- length(bin_labels)
  bin_results <- vector("list", B)
  
  for (b_idx in seq_len(B)) {
    this_bin <- bin_labels[b_idx]
    sub_df   <- df[bin_id == this_bin]
    
    if (nrow(sub_df)==0) {
      bin_results[[b_idx]] <- data.frame(
        bin       = this_bin,
        bin_range = paste0("[",bin_breaks[b_idx],",",bin_breaks[b_idx+1],")"),
        n_snps    = 0,
        n_blocks  = 0,
        mean_eta  = NA_real_,
        se        = NA_real_,
        lower_ci  = NA_real_,
        upper_ci  = NA_real_
      )
      next
    }
    # Filter blocks
    counts_block <- sub_df[,.N, by=block_id]
    valid_blocks <- counts_block[N>=min_snps]$block_id
    sub_df <- sub_df[block_id %in% valid_blocks]
    
    if (nrow(sub_df)==0) {
      bin_results[[b_idx]] <- data.frame(
        bin       = this_bin,
        bin_range = paste0("[",bin_breaks[b_idx],",",bin_breaks[b_idx+1],")"),
        n_snps    = 0,
        n_blocks  = 0,
        mean_eta  = NA_real_,
        se        = NA_real_,
        lower_ci  = NA_real_,
        upper_ci  = NA_real_
      )
      next
    }
    
    blocks_in_bin <- unique(sub_df$block_id)
    n_snps_bin    <- nrow(sub_df)
    n_blocks_bin  <- length(blocks_in_bin)
    
    # pick lowest pval
    block_map <- split(sub_df, sub_df$block_id)
    chosen_etas <- sapply(block_map, function(bd){
      best_row <- bd[which.min(bd$pval)]
      best_row$eta
    })
    point_est <- sum(chosen_etas)/sum(abs(chosen_etas))
    
    # block bootstrap
    block_ids <- blocks_in_bin
    K <- length(block_ids)
    boot_vals <- numeric(n_boot)
    for (i in seq_len(n_boot)) {
      bs_blocks <- sample(block_ids, size=K, replace=TRUE)
      chosen_bs <- c()
      for (bk in bs_blocks) {
        block_sub <- block_map[[as.character(bk)]]
        if (!is.null(block_sub)) {
          best_row <- block_sub[which.min(block_sub$pval)]
          chosen_bs <- c(chosen_bs, best_row$eta)
        }
      }
      if (length(chosen_bs)==0) {
        boot_vals[i] <- NA_real_
      } else {
        boot_vals[i] <- sum(chosen_bs)/sum(abs(chosen_bs))
      }
    }
    sd_boot <- sd(boot_vals, na.rm=TRUE)
    l_ci <- point_est - 1.96*sd_boot
    u_ci <- point_est + 1.96*sd_boot
    
    bin_results[[b_idx]] <- data.frame(
      bin       = this_bin,
      bin_range = paste0("[",bin_breaks[b_idx],",",bin_breaks[b_idx+1],")"),
      n_snps    = n_snps_bin,
      n_blocks  = n_blocks_bin,
      mean_eta  = point_est,
      se        = sd_boot,
      lower_ci  = l_ci,
      upper_ci  = u_ci
    )
  }
  
  bin_df <- do.call(rbind, bin_results)
  rownames(bin_df) <- NULL
  return(bin_df)
}

###############################################################################
# 3) Functions for THRESHOLD mode
###############################################################################

calc_threshold_random <- function(df, freq_type="MAF", T_low=0.01, T_high=NA, N=1000, min_snps=10) {
  # Filter
  if (freq_type=="MAF") {
    sub_df <- df[af < T_low]
  } else {
    # DAF
    if (!is.na(T_high)) {
      sub_df <- df[af < T_low | af > T_high]
    } else {
      sub_df <- df[af < T_low]
    }
  }
  if (nrow(sub_df)==0) {
    return(data.frame(
      mean_eta  = NA_real_, se=NA_real_, lower_ci=NA_real_, upper_ci=NA_real_,
      n_snps=0, n_blocks=0
    ))
  }
  
  # filter blocks
  counts_block <- sub_df[,.N, by=block_id]
  valid_blocks <- counts_block[N>=min_snps]$block_id
  sub_df <- sub_df[block_id %in% valid_blocks]
  
  if (nrow(sub_df)==0) {
    return(data.frame(
      mean_eta=NA_real_, se=NA_real_, lower_ci=NA_real_, upper_ci=NA_real_,
      n_snps=0, n_blocks=0
    ))
  }
  
  blocks_in_data <- unique(sub_df$block_id)
  n_snps   <- nrow(sub_df)
  n_blocks <- length(blocks_in_data)
  block_map <- split(sub_df$eta, sub_df$block_id)
  
  replicate_etas <- numeric(N)
  for (i in seq_len(N)) {
    chosen <- sapply(block_map, function(etas) sample(etas,1))
    replicate_etas[i] <- sum(chosen)/sum(abs(chosen))
  }
  mean_eta <- mean(replicate_etas)
  sd_eta   <- sd(replicate_etas)
  l_ci     <- mean_eta - 1.96*sd_eta
  u_ci     <- mean_eta + 1.96*sd_eta
  
  out <- data.frame(
    mean_eta  = mean_eta,
    se        = sd_eta,
    lower_ci  = l_ci,
    upper_ci  = u_ci,
    n_snps    = n_snps,
    n_blocks  = n_blocks
  )
  return(out)
}

calc_threshold_sig <- function(df, freq_type="MAF", T_low=0.01, T_high=NA, n_boot=1000, min_snps=10) {
  if (freq_type=="MAF") {
    sub_df <- df[af < T_low]
  } else {
    if (!is.na(T_high)) {
      sub_df <- df[af < T_low | af > T_high]
    } else {
      sub_df <- df[af < T_low]
    }
  }
  if (nrow(sub_df)==0) {
    return(data.frame(
      mean_eta=NA_real_, se=NA_real_, lower_ci=NA_real_, upper_ci=NA_real_,
      n_snps=0, n_blocks=0
    ))
  }
  
  counts_block <- sub_df[,.N, by=block_id]
  valid_blocks <- counts_block[N>=min_snps]$block_id
  sub_df <- sub_df[block_id %in% valid_blocks]
  
  if (nrow(sub_df)==0) {
    return(data.frame(
      mean_eta=NA_real_, se=NA_real_, lower_ci=NA_real_, upper_ci=NA_real_,
      n_snps=0, n_blocks=0
    ))
  }
  
  blocks_in_data <- unique(sub_df$block_id)
  n_snps   <- nrow(sub_df)
  n_blocks <- length(blocks_in_data)
  block_map <- split(sub_df, sub_df$block_id)
  
  # pick lowest p
  chosen_etas <- sapply(block_map, function(bd){
    best_row <- bd[which.min(bd$pval)]
    best_row$eta
  })
  point_est <- sum(chosen_etas)/sum(abs(chosen_etas))
  boot_vals <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    bs_blocks <- sample(blocks_in_data, size=n_blocks, replace=TRUE)
    chosen_bs <- c()
    for (bk in bs_blocks) {
      block_sub <- block_map[[as.character(bk)]]
      if (!is.null(block_sub)) {
        best_row <- block_sub[which.min(block_sub$pval)]
        chosen_bs <- c(chosen_bs, best_row$eta)
      }
    }
    if (length(chosen_bs)==0) {
      boot_vals[i] <- NA_real_
    } else {
      boot_vals[i] <- sum(chosen_bs)/sum(abs(chosen_bs))
    }
  }
  sd_boot <- sd(boot_vals, na.rm=TRUE)
  l_ci <- point_est - 1.96*sd_boot
  u_ci <- point_est + 1.96*sd_boot
  
  out <- data.frame(
    mean_eta  = point_est,
    se        = sd_boot,
    lower_ci  = l_ci,
    upper_ci  = u_ci,
    n_snps    = n_snps,
    n_blocks  = n_blocks
  )
  return(out)
}

###############################################################################
# 4) Main: Read data and run the chosen approach
###############################################################################

cat("Reading file...\n")

# We'll pick columns to read depending on type:
if (type_arg=="MAF") {
  keep_cols <- c("block_id","MAF","eta","pval")
} else if (type_arg=="DAF") {
  keep_cols <- c("block_id","DAF","eta","pval")
} else {
  stop("type_arg must be MAF or DAF.")
}

df <- fread(file_arg, select=keep_cols)

# rename frequency col -> "af"
if (type_arg=="MAF") {
  if (!"MAF" %in% names(df)) {
    stop("No 'MAF' column found after renaming.")
  }
  df[, af := MAF]
} else {
  if (!"DAF" %in% names(df)) {
    stop("No 'DAF' column found after renaming.")
  }
  df[, af := DAF]
}

cat("Running method...\n")

if (mode_arg=="bin") {
  # BIN approach
  if (method_arg=="random") {
    if (type_arg=="MAF") {
      bin_df <- calc_bin_stats_random(df, bin_breaks_maf, N=1000, min_snps=10)
    } else {
      bin_df <- calc_bin_stats_random(df, bin_breaks_daf, N=1000, min_snps=10)
    }
    out_binfile <- paste0("bin_", trait_name, "_", type_arg, "_random.txt")
    fwrite(
      bin_df %>% mutate(trait=trait_name, type=type_arg, method="random"),
      file=out_binfile, sep="\t", quote=FALSE, row.names=FALSE
    )
    cat("Wrote bin-level file:", out_binfile, "\n")
    
  } else if (method_arg=="sig") {
    if (type_arg=="MAF") {
      bin_df <- calc_bin_stats_sig(df, bin_breaks_maf, n_boot=1000, min_snps=10)
    } else {
      bin_df <- calc_bin_stats_sig(df, bin_breaks_daf, n_boot=1000, min_snps=10)
    }
    out_binfile <- paste0("bin_", trait_name, "_", type_arg, "_sig.txt")
    fwrite(
      bin_df %>% mutate(trait=trait_name, type=type_arg, method="sig"),
      file=out_binfile, sep="\t", quote=FALSE, row.names=FALSE
    )
    cat("Wrote bin-level file:", out_binfile, "\n")
    
  } else {
    stop("method_arg must be random or sig for bin mode.")
  }
  
} else if (mode_arg=="threshold") {
  # THRESHOLD approach
  if (is.na(T_low) & is.na(T_high)) {
    stop("For threshold mode, please specify at least --low or --high.")
  }
  
  # We'll produce a single row of results with (mean_eta, se, lower_ci, upper_ci, n_snps, n_blocks)
  # We'll incorporate T_low / T_high in the filename
  low_str  <- if(!is.na(T_low)) paste0("_",T_low) else ""
  high_str <- if(!is.na(T_high)) paste0("_",T_high) else ""
  
  if (method_arg=="random") {
    out_df <- calc_threshold_random(
      df,
      freq_type=type_arg,
      T_low=if(is.na(T_low)) 0.0 else T_low,
      T_high=T_high,
      N=1000,
      min_snps=10
    )
    out_df$trait   <- trait_name
    out_df$type    <- type_arg
    out_df$method  <- "random"
    out_df$T_low   <- T_low
    out_df$T_high  <- T_high
    
    out_file <- paste0("threshold_", trait_name, "_", type_arg, "_random", low_str, high_str, ".txt")
    fwrite(out_df, file=out_file, sep="\t", quote=FALSE, row.names=FALSE)
    cat("Wrote threshold summary file:", out_file, "\n")
    
  } else if (method_arg=="sig") {
    out_df <- calc_threshold_sig(
      df,
      freq_type=type_arg,
      T_low=if(is.na(T_low)) 0.0 else T_low,
      T_high=T_high,
      n_boot=1000,
      min_snps=10
    )
    out_df$trait   <- trait_name
    out_df$type    <- type_arg
    out_df$method  <- "sig"
    out_df$T_low   <- T_low
    out_df$T_high  <- T_high
    
    out_file <- paste0("threshold_", trait_name, "_", type_arg, "_sig", low_str, high_str, ".txt")
    fwrite(out_df, file=out_file, sep="\t", quote=FALSE, row.names=FALSE)
    cat("Wrote threshold summary file:", out_file, "\n")
    
  } else {
    stop("method_arg must be random or sig for threshold mode.")
  }
  
} else {
  stop("mode_arg must be bin or threshold.")
}

cat("Done.\n")

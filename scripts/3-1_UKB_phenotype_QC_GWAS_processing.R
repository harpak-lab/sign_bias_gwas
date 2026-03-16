##########################################################
## Sample QC / Covariates / Prepare phenotypes for GWAS ##
##########################################################

#libraries/wd
library(dplyr)
library(data.table)
setwd("ukb_group_gwas/phenos/")

ukb      <- fread("UKB.phenos.all.csv")
ukb_qc   <- fread("UKB.qc.csv")
ukb      <- merge(ukb, ukb_qc, by = "eid")

gen_samples <- fread("ukb61666_imp_chr1_v3_s487280.sample", 
  sep = " ")
gen_samples <- gen_samples[-1, ]

ukb <- ukb[eid %in% gen_samples$ID_1]
setDT(ukb)


######################## Sample QC #######################

##Coerce numeric QC fields 
qc_numeric <- c(
  "31-0.0",      # reported sex
  "21000-0.0",   # self-reported ethnicity (instance 0)
  "21022-0.0",   # age at assessment
  "22000-0.0",   # batch
  "22001-0.0",   # genetic sex
  "22003-0.0",   # heterozygosity
  "22004-0.0",   # heterozygosity, PCA corrected
  "22005-0.0",   # missingness
  "22006-0.0",   # genetic ethnic grouping
  "22019-0.0",   # sex chromosome aneuploidy
  "22020-0.0",   # used in genetic PCs
  "22021-0.0",   # kinship summary
  "22022-0.0",   # X intensity
  "22023-0.0",   # Y intensity
  "22024-0.0",   # DNA concentration
  "22025-0.0",   # Cluster.CR
  "22026-0.0",   # dQC
  "22027-0.0",   # outliers for het/missing
  "22028-0.0",   # use in phasing chr1-22
  "22029-0.0",   # use in phasing chrX
  "22030-0.0"    # use in phasing chrXY
)

qc_numeric <- intersect(qc_numeric, names(ukb))
ukb[, (qc_numeric) := lapply(.SD, as.numeric), .SDcols = qc_numeric]

## QC flags

# Heterozygosity/missing outliers: 22027 
ukb[, het_missing_outliers :=
      fifelse(`22027-0.0` == 1, TRUE, FALSE, na = FALSE)]

# Sex chromosome aneuploidy: 22019 == 1
ukb[, putative_sex_chromosome_aneuploidy :=
      fifelse(`22019-0.0` == 1, TRUE, FALSE, na = FALSE)]

# Kinship flags from 22021
# -1: excluded from kinship inference
#  0: no kinship
#  1: >= 1 relative
# 10: >= 10 relatives

ukb[, excluded_from_kinship_inference :=
      fifelse(`22021-0.0` == -1, TRUE, FALSE, na = FALSE)]

ukb[, excess_relatives :=
      fifelse(`22021-0.0` == 10, TRUE, FALSE, na = FALSE)]

ukb[, in_kinship_table :=
      fifelse(!is.na(`22021-0.0`) & `22021-0.0` != 0, TRUE, FALSE)]

# Used in PCA calculation: 22020 == 1
ukb[, used_in_pca_calculation :=
      fifelse(`22020-0.0` == 1, TRUE, FALSE, na = FALSE)]

# Phasing flags
ukb[, in_Phasing_Input_chr1_22 :=
      fifelse(`22028-0.0` == 1, TRUE, FALSE, na = FALSE)]
ukb[, in_Phasing_Input_chrX :=
      fifelse(`22029-0.0` == 1, TRUE, FALSE, na = FALSE)]
ukb[, in_Phasing_Input_chrXY :=
      fifelse(`22030-0.0` == 1, TRUE, FALSE, na = FALSE)]

# Genetic White British ancestry subset (UKB genetic group 1)
ukb[, in_white_British_ancestry_subset :=
      fifelse(`22006-0.0` == 1, TRUE, FALSE, na = FALSE)]

## Age/ancestry flags

# Age at assessment
ukb[, age := as.numeric(`21022-0.0`)]
ukb[, is_in_phenotype_data := !is.na(age)]

# Self-reported ethnicity (instance 0)
ukb[, ethnic := `21000-0.0`]

# Self-reported "white" 
# 1    = White
# 1001 = White British
# 1002 = White Irish
# 1003 = Any other white background
ukb[, white       := ethnic %in% c(1, 1001, 1002, 1003)]
ukb[, ethnic_miss := ethnic %in% c(-3, -1) | is.na(ethnic)]

## 4. Rename PCs 22009-0.x → PC1..PC40
pc_cols <- grep("^22009-0\\.", names(ukb), value = TRUE)

setnames(
  ukb,
  old = pc_cols,
  new = paste0("PC", seq_along(pc_cols))
)

pc_nams <- paste0("PC", seq_along(pc_cols))   # PC1..PC40

## 5. Base QC mask 
ukb[, base_qc_pass := (
  !het_missing_outliers &
    !putative_sex_chromosome_aneuploidy &
    !excluded_from_kinship_inference &
    !excess_relatives &
    used_in_pca_calculation
)]

##British-ancestry ellipse

# Use only individuals in the curated WB subset *and* passing base QC
wb_idx <- which(ukb$in_white_British_ancestry_subset & ukb$base_qc_pass)

ells     <- 6      # first 6 PCs
sds_brit <- 7      # 7 SD radius as per Neale
pc_use   <- paste0("PC", 1:ells)

mm_brit <- sapply(pc_use, function(p) mean(ukb[wb_idx, get(p)], na.rm = TRUE))
ss_brit <- sapply(pc_use, function(p) sd(ukb[wb_idx, get(p)],   na.rm = TRUE))

dd_brit <- rep(0, nrow(ukb))
for (i in seq_len(ells)) {
  pc_vals <- ukb[[pc_use[i]]]
  dd_brit <- dd_brit + (pc_vals - mm_brit[i])^2 / (ss_brit[i]^2)
}
ukb[, dd_brit := dd_brit]

# Neale's EUR selection:
# - within 7 SD ellipse around WB cluster
ukb[, eur_select := (dd_brit < sds_brit^2) & (white | ethnic_miss)]
ukb[, is_in_ancestry_subset := eur_select]

## Final QC sample
ukb[, final_gwas_sample := (
  base_qc_pass &
    in_Phasing_Input_chr1_22 &
    in_Phasing_Input_chrX &
    is_in_phenotype_data &
    is_in_ancestry_subset &
    used_in_pca_calculation
)]

ukb_final <- ukb[final_gwas_sample == TRUE]

## output samples (keep file)
keep_qc <- ukb_final[, .(
  FID = eid,
  IID = eid
)]

fwrite(
  keep_qc,
  file      = "ukb_qc_samples.txt",
  sep       = "\t",
  col.names = TRUE,
  quote     = FALSE
)

######################## Covariates #######################
# must have "age"
if (!"age" %in% names(ukb_final)) {
  ukb_final[, age := as.numeric(`21022-0.0`)]
}

# sex (numeric)
ukb_final[, sex := as.numeric(`22001-0.0`)]

# Drop any samples with missing sex/age 
ukb_final <- ukb_final[!is.na(sex) & !is.na(age)]

# Age-derived terms
ukb_final[, age2      := age^2]
ukb_final[, sex_age   := sex * age]
ukb_final[, sex_age2  := sex * age2]

pc_use <- paste0("PC", 1:20)

# Add FID/IID columns
ukb_final[, FID := eid]
ukb_final[, IID := eid]

covar_cols <- c(
  "FID", "IID",
  "sex",
  "age", "age2",
  "sex_age", "sex_age2",
  pc_use
)

covar_dt <- ukb_final[, ..covar_cols]

fwrite(
  covar_dt,
  file      = "ukb.covariates.txt",
  sep       = "\t",
  col.names = TRUE,
  quote     = FALSE
)

################## Build binary/quant phenotypes #################

ukb_binary <- as.data.table(ukb_final)

## Identify main ICD10 columns (41202-*) 
icd_main_cols <- grep("^41202-", names(ukb_binary), value = TRUE)

icd_main_long <- melt(
  ukb_binary,
  id.vars       = "eid",
  measure.vars  = icd_main_cols,
  variable.name = "source",
  value.name    = "icd10_raw",
  variable.factor = FALSE
)

# Keep only non-missing, non-empty codes; normalize case/whitespace
icd_main_long <- icd_main_long[!is.na(icd10_raw) & icd10_raw != ""]
icd_main_long[, icd10_raw := toupper(trimws(as.character(icd10_raw)))]

# 3-character ICD root, e.g. "G30", "E10", "E11", "F20"
icd_main_long[, icd3 := substr(icd10_raw, 1, 3)]

## 3. Define case sets for each endpoint 

# Alzheimer's disease (G6_ALZHEIMER): any G30* in main ICD10
ad_ids <- icd_main_long[icd3 == "G30", unique(eid)]

# Asthma (J45_ASTHMA): any G30* in main ICD10
ast_ids <- icd_main_long[icd3 == "J45", unique(eid)]

# Type 1 diabetes (E4_DM1): any E10* in main ICD10
t1d_ids <- icd_main_long[icd3 == "E10", unique(eid)]

# Type 2 diabetes (E4_DM2): any E11* in main ICD10
t2d_ids_raw <- icd_main_long[icd3 == "E11", unique(eid)]

# Hierarchical rule: T1D has priority; remove T1D from T2D
t2d_ids <- setdiff(t2d_ids_raw, t1d_ids)

# Schizophrenia (SCZ): main ICD10 F20* (Neale: Diagnoses – main ICD10: F20)
scz_ids <- icd_main_long[icd3 == "F20", unique(eid)]

## Attach binary phenotypes 
ukb_binary[, AD  := as.integer(eid %in% ad_ids)]
ukb_binary[, ASTHMA  := as.integer(eid %in% ast_ids)]
ukb_binary[, T1D := as.integer(eid %in% t1d_ids)]
ukb_binary[, T2D := as.integer(eid %in% t2d_ids)]
ukb_binary[, SCZ := as.integer(eid %in% scz_ids)]

##  counts
summary_counts <- rbindlist(list(
  ukb_binary[, .(trait = "AD",
         n_total    = .N,
         n_cases    = sum(AD,  na.rm = TRUE),
         n_controls = sum(AD  == 0, na.rm = TRUE))],
  ukb_binary[, .(trait = "T1D",
         n_total    = .N,
         n_cases    = sum(T1D, na.rm = TRUE),
         n_controls = sum(T1D == 0, na.rm = TRUE))],
  ukb_binary[, .(trait = "T2D",
         n_total    = .N,
         n_cases    = sum(T2D, na.rm = TRUE),
         n_controls = sum(T2D == 0, na.rm = TRUE))],
  ukb_binary[, .(trait = "SCZ",
         n_total    = .N,
         n_cases    = sum(SCZ, na.rm = TRUE),
         n_controls = sum(SCZ == 0, na.rm = TRUE))],
  ukb_binary[, .(trait = "ASTHMA",
                 n_total    = .N,
                 n_cases    = sum(ASTHMA, na.rm = TRUE),
                 n_controls = sum(ASTHMA == 0, na.rm = TRUE))]
), use.names = TRUE)

summary_counts$neale_cases <- c(119,583,888,
  198,1693)
summary_counts$neale_controls <- c(361075,360611,360306,
  360996,359501)


print(summary_counts)

## Update pheno table
ukb_final <- ukb_binary

##integrate quant traits
ukb_quant <- as.data.table(ukb_final)

# Map quant traits
quant_map <- c(
  height    = "50-0.0",     # standing height
  bmi       = "21001-0.0",  # BMI
  weight    = "21002-0.0",  # weight
  mono_pct  = "30190-0.0",  # monocyte %
  neutro_pct= "30200-0.0",  # neutrophil %
  baso_pct  = "30220-0.0",   # basophil %
  wbc = "30000-0.0", # white blood cell count
  rbc = "30010-0.0", # red blood cell count
  mch = "30050-0.0" # mean corpuscular hemoglobin
)

for (new_name in names(quant_map)) {
  old_col <- quant_map[[new_name]]
  if (!old_col %in% names(ukb_quant)) {
    warning(sprintf("Column %s not found in ukb_final; %s will be NA.", old_col, new_name))
    ukb_quant[, (new_name) := NA_real_]
  } else {
    ukb_quant[, (new_name) := as.numeric(get(old_col))]
  }
}

# Update pheno table
ukb_final <- ukb_quant

# Export phenotypes (PLINK)

# add +2 to each binary value
binary_traits <- c("AD", "T1D", 
                   "T2D", "SCZ", 
                   "ASTHMA")  

ukb_final[, (binary_traits) := lapply(.SD, function(x) {
  x <- as.integer(x)                # ensure numeric/integer
  x[!is.na(x)] <- x[!is.na(x)] + 2L # add 2 only where not NA
  x
}), .SDcols = binary_traits]

#define all phenos
pheno_cols <- c(
  "AD", "T1D", "T2D", "SCZ", "ASTHMA",       # binary
  "height", "bmi", "weight",                # quantitative
  "mono_pct", "neutro_pct", "baso_pct",     # blood percentages
  "wbc", "rbc", "mch"
)

# Construct PLINK table
pheno_plink <- ukb_final[
  , c("IID", pheno_cols),
  with = FALSE
]

# FID = IID = eid
setnames(pheno_plink, "IID", "eid")   # temporary rename
pheno_plink[, IID := eid]
pheno_plink[, FID := eid]

setcolorder(pheno_plink, c("FID", "IID", pheno_cols))

# Drop eid and export
pheno_plink[, eid := NULL]

fwrite(
  pheno_plink,
  file = "ukb_phenotypes.txt",
  sep = "\t",
  quote = FALSE,
  na = "NA")

################## Split QC samples #################
qc<- keep_qc

#add a random number and sort
qc[, rand := runif(.N)]
setorder(qc, rand)
qc[, rand := NULL]

#assign first half to group 1, second half to group 2
n <- nrow(qc)
qc[, split2 := ifelse(1:.N <= n/2, 1L, 2L)]

table(qc$split2)
# should be roughly n/2, n/2

# write out two split files, preserving original columns
orig_cols <- setdiff(names(qc), "split2")

qc_split1 <- qc[split2 == 1L, ..orig_cols]
qc_split2 <- qc[split2 == 2L, ..orig_cols]

fwrite(qc_split1, "ukb_samples_group1.txt",
       sep = "\t", col.names = TRUE, quote = FALSE)
fwrite(qc_split2, "ukb_samples_group2.txt",
       sep = "\t", col.names = TRUE, quote = FALSE)

### Split into 3
qc <- keep_qc

# Add a random number and sort
qc[, rand := runif(.N)]
setorder(qc, rand)
qc[, rand := NULL]

# Assign rows to 3 groups: 1, 2, 3 (roughly equal sizes)
n <- nrow(qc)
qc[, split3 := ((1:.N - 1L) %% 3L) + 1L]

table(qc$split3)
# should be roughly n/3 in each of 1, 2, 3

# Write out three split files, preserving original columns
orig_cols <- setdiff(names(qc), "split3")

qc_split1 <- qc[split3 == 1L, ..orig_cols]
qc_split2 <- qc[split3 == 2L, ..orig_cols]
qc_split3 <- qc[split3 == 3L, ..orig_cols]

fwrite(qc_split1, "ukb_samples_smallgroup1.txt",
       sep = "\t", col.names = TRUE, quote = FALSE)
fwrite(qc_split2, "ukb_samples_smallgroup2.txt",
       sep = "\t", col.names = TRUE, quote = FALSE)
fwrite(qc_split3, "ukb_samples_smallgroup3.txt",
       sep = "\t", col.names = TRUE, quote = FALSE)

################## 

# Binary traits in ukb_final (already recoded: 0->2, 1->3)
binary_traits <- c("AD", "T1D", "T2D", "SCZ", "ASTHMA")

# Helper: summarize one split
summarize_split <- function(keep_file, split_label) {
  # Read keep file: FID / IID (IID = eid)
  ids <- fread(keep_file)
  
  # Make an 'eid' column to join on
  if (!"eid" %in% names(ids)) {
    ids[, eid := IID]
  }
  
  # Restrict ukb_final to this split
  split_dt <- ukb_final[eid %in% ids$eid]
  
  # Build summary for each binary trait
  rbindlist(lapply(binary_traits, function(tr) {
    v <- split_dt[[tr]]
    data.table(
      split      = split_label,
      trait      = tr,
      n_total    = length(v),
      n_cases    = sum(v == 3L, na.rm = TRUE),  # 3 = case (1+2)
      n_controls = sum(v == 2L, na.rm = TRUE),  # 2 = control (0+2)
      n_missing  = sum(is.na(v))
    )
  }))
}

# Collect summaries for all 5 splits
summary_splits <- rbindlist(list(
  summarize_split("ukb_samples_group1.txt",       "big_group1"),
  summarize_split("ukb_samples_group2.txt",       "big_group2"),
  summarize_split("ukb_samples_smallgroup1.txt",  "small_group1"),
  summarize_split("ukb_samples_smallgroup2.txt",  "small_group2"),
  summarize_split("ukb_samples_smallgroup3.txt",  "small_group3")
), use.names = TRUE)

print(summary_splits)
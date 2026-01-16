############################################################
## AoU phenotype/sample QC  
############################################################

library(data.table)

bucket <- Sys.getenv("WORKSPACE_BUCKET")
stopifnot(nzchar(bucket))

PHENO_DIR <- "phenos"
dir.create(PHENO_DIR, showWarnings = FALSE, recursive = TRUE)

gs_cp <- function(src, dest = ".") {
  system2("gsutil", c("-m", "cp", src, dest))
}

trim_id <- function(x) trimws(as.character(x))

calc_central_moments <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0L) return(c(mu0=NA_real_, mu1=NA_real_, mu2=NA_real_, mu3=NA_real_, mu3_std=NA_real_))
  m  <- mean(x)
  xc <- x - m
  mu2 <- mean(xc^2)
  mu3 <- mean(xc^3)
  skew <- if (is.finite(mu2) && mu2 > 0) mu3 / (mu2^(3/2)) else NA_real_
  c(mu0=1, mu1=0, mu2=mu2, mu3=mu3, mu3_std=skew)
}

# Writes person_id/value phenotype
write_pheno <- function(DT, id_col, val_col, keep_ids, out) {
  tmp <- DT[, .(person_id = get(id_col), value = get(val_col))]
  tmp[, person_id := trim_id(person_id)]
  tmp <- tmp[person_id %chin% keep_ids & !is.na(value)]
  fwrite(tmp, out, sep = "\t")
  invisible(tmp)
}

# Writes PLINK phenotype: FID IID PHENO
write_pheno_plink <- function(pheno_dt, out_plink) {
  stopifnot(all(c("person_id","value") %in% names(pheno_dt)))
  pl <- pheno_dt[, .(FID = person_id, IID = person_id, PHENO = value)]
  fwrite(pl, out_plink, sep = "\t", col.names = TRUE)
  invisible(pl)
}

# Mean-center covariates 
center_covariates <- function(cov_dt, id_col = "person_id") {
  covc <- copy(cov_dt)
  num_cols <- names(covc)[vapply(covc, is.numeric, logical(1L))]
  num_cols <- setdiff(num_cols, id_col)
  if (length(num_cols) == 0L) return(covc)

  means <- covc[, lapply(.SD, mean, na.rm = TRUE), .SDcols = num_cols]
  for (j in num_cols) {
    covc[, (j) := get(j) - means[[j]]]
  }
  covc
}


## Download inputs
remote_files <- paste0(bucket, "/", c(
  "aou_gwas/phecodeX_files/mcc2_phecodex_table_v8.csv",
  "aou_gwas/ancestry_table.csv",
  "aou_gwas/demographics_table.csv",
  "aou_gwas/participant_PCs.csv",
  "aou_gwas/samples_relatedness_flagged_samples.tsv",
  "aou_gwas/aou.v8.classified-rf.wb-classify.report-any.fam",
  "aou_gwas/pheno/lab_measures_wbc.tsv",
  "aou_gwas/pheno/physical_measures_hwb.tsv"
))
gs_cp(remote_files, ".")

## Build binary phenos subset
bin_cols <- c("person_id", "sex", "RE_475","EM_202.1",
  "EM_202.2","MB_287.1","NS_328.11")
binary_table <- fread("mcc2_phecodex_table_v8.csv", select = bin_cols)
binary_table[, person_id := trim_id(person_id)]
fwrite(binary_table, file.path(PHENO_DIR, "binary_phenos.csv"))


## Read other inputs
demo_table <- fread("demographics_table.csv", select = 1:3)
setnames(demo_table, c("person_id","age","sex"))
demo_table[, person_id := trim_id(person_id)]

pcs <- fread("participant_PCs.csv", select = c(1, 5:20))
pcs[, person_id := trim_id(person_id)]

relatedness <- fread("samples_relatedness_flagged_samples.tsv", header = FALSE)
setnames(relatedness, "person_id")
relatedness[, person_id := trim_id(person_id)]

wb_fam <- fread("aou.v8.classified-rf.wb-classify.report-any.fam", header = FALSE)
wb_ukb_samples <- wb_fam[, .(person_id = trim_id(V2))]

physical <- fread("physical_measures_hwb.tsv", sep = "\t")
physical[, person_id := trim_id(person_id)]

lab_wb <- fread("lab_measures_wbc.tsv", sep = "\t")
lab_wb[, person_id := trim_id(person_id)]


## Sample QC: unrelated IDs
related_ids <- unique(relatedness[!is.na(person_id), person_id])
wb_ids      <- unique(wb_ukb_samples[!is.na(person_id), person_id])
unrelated_ids <- setdiff(wb_ids, related_ids)

# keep file
wb_ukb_plink <- data.table(`#IID` = unrelated_ids)
fwrite(wb_ukb_plink, "wb_samples_filtered.txt", sep = "\t")


## Covariates (centered)
covariates <- merge(demo_table, pcs, by = "person_id", all = FALSE)
covariates <- covariates[sex %chin% c("Male","Female") & person_id %chin% unrelated_ids]

covariates[, age2 := age^2]
covariates[, sex := fifelse(sex == "Male", 1L, 0L)]
covariates[, `:=`(
  sex_age  = age  * sex,
  sex_age2 = age2 * sex
)]

cov_centered <- center_covariates(covariates, id_col = "person_id")
fwrite(cov_centered, file.path(PHENO_DIR, "aou_covariates.centered.tsv"), sep = "\t")



## Pheno QC: write to phenos
phys_specs <- list(
  height = list(val = "height_cm", out = file.path(PHENO_DIR, "aou_height.tsv")),
  weight = list(val = "weight_kg", out = file.path(PHENO_DIR, "aou_weight.tsv")),
  bmi    = list(val = "bmi",       out = file.path(PHENO_DIR, "aou_bmi.tsv"))
)

phys_out <- list()
for (nm in names(phys_specs)) {
  s <- phys_specs[[nm]]
  phys_out[[nm]] <- write_pheno(physical, "person_id", s$val, unrelated_ids, s$out)
  write_pheno_plink(phys_out[[nm]], sub("\\.tsv$", ".plink.tsv", s$out))
}

# Lab
lab_map <- c(
  monocyte_percentage    = file.path(PHENO_DIR, "aou_monocyte_percentage.tsv"),
  basophil_percentage    = file.path(PHENO_DIR, "aou_basophil_percentage.tsv"),
  neutrophil_percentage  = file.path(PHENO_DIR, "aou_neutrophil_percentage.tsv")
)

lab_out <- list()
for (lab_name in names(lab_map)) {
  tmp <- lab_wb[lab == lab_name, .(person_id, median)]
  setnames(tmp, "median", "value")
  tmp[, person_id := trim_id(person_id)]
  tmp <- tmp[person_id %chin% unrelated_ids & !is.na(value)]
  fwrite(tmp, lab_map[[lab_name]], sep = "\t")
  write_pheno_plink(tmp, sub("\\.tsv$", ".plink.tsv", lab_map[[lab_name]]))
  lab_out[[lab_name]] <- tmp
}

# Binary
bin_specs <- list(
  asthma        = list(col = "RE_475",    out = file.path(PHENO_DIR, "aou_asthma.tsv")),
  t1d           = list(col = "EM_202.1",  out = file.path(PHENO_DIR, "aou_t1d.tsv")),
  t2d           = list(col = "EM_202.2",  out = file.path(PHENO_DIR, "aou_t2d.tsv")),
  schizophrenia = list(col = "MB_287.1",  out = file.path(PHENO_DIR, "aou_schizophrenia.tsv")),
  alzheimers    = list(col = "NS_328.11", out = file.path(PHENO_DIR, "aou_alzheimers.tsv"))
)

bin_out <- list()
for (nm in names(bin_specs)) {
  s <- bin_specs[[nm]]
  tmp <- binary_table[, .(person_id, value = get(s$col))]
  tmp <- tmp[person_id %chin% unrelated_ids & !is.na(value)]
  fwrite(tmp, s$out, sep = "\t")
  write_pheno_plink(tmp, sub("\\.tsv$", ".plink.tsv", s$out))
  bin_out[[nm]] <- tmp
}


## Moments (write to phenos)
trait_vectors <- list(
  "Height"                 = phys_out$height$value,
  "Weight"                 = phys_out$weight$value,
  "BMI"                    = phys_out$bmi$value,
  "Monocyte percentage"    = lab_out$monocyte_percentage$value,
  "Basophil percentage"    = lab_out$basophil_percentage$value,
  "Neutrophil percentage"  = lab_out$neutrophil_percentage$value,
  "Asthma"                 = bin_out$asthma$value,
  "Type 1 diabetes"        = bin_out$t1d$value,
  "Type 2 diabetes"        = bin_out$t2d$value,
  "Schizophrenia"          = bin_out$schizophrenia$value,
  "Alzheimer's disease"    = bin_out$alzheimers$value
)

aou_moments <- rbindlist(lapply(names(trait_vectors), function(tr) {
  moms <- calc_central_moments(trait_vectors[[tr]])
  data.table(Trait = tr,
             mu0 = moms["mu0"], mu1 = moms["mu1"], mu2 = moms["mu2"], mu3 = moms["mu3"], mu3_std = moms["mu3_std"])
}), use.names = TRUE)

fwrite(aou_moments, file.path(PHENO_DIR, "aou_moments.tsv"), sep = "\t")

#Upload outputs to bucket
out_bucket_map <- list(
  file.path(PHENO_DIR, "aou_covariates.centered.tsv") = paste0(bucket, "/aou_gwas/"),
  "wb_samples_filtered.txt"                           = paste0(bucket, "/aou_gwas/"),

  # phenos + plink phenos
  file.path(PHENO_DIR, "aou_height.tsv")              = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_height.plink.tsv")        = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_weight.tsv")              = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_weight.plink.tsv")        = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_bmi.tsv")                 = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_bmi.plink.tsv")           = paste0(bucket, "/aou_gwas/pheno/"),

  file.path(PHENO_DIR, "aou_monocyte_percentage.tsv")       = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_monocyte_percentage.plink.tsv") = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_basophil_percentage.tsv")       = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_basophil_percentage.plink.tsv") = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_neutrophil_percentage.tsv")     = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_neutrophil_percentage.plink.tsv") = paste0(bucket, "/aou_gwas/pheno/"),

  file.path(PHENO_DIR, "aou_asthma.tsv")              = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_asthma.plink.tsv")        = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_t1d.tsv")                 = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_t1d.plink.tsv")           = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_t2d.tsv")                 = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_t2d.plink.tsv")           = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_schizophrenia.tsv")       = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_schizophrenia.plink.tsv") = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_alzheimers.tsv")          = paste0(bucket, "/aou_gwas/pheno/"),
  file.path(PHENO_DIR, "aou_alzheimers.plink.tsv")    = paste0(bucket, "/aou_gwas/pheno/"),

  file.path(PHENO_DIR, "aou_moments.tsv")             = paste0(bucket, "/aou_gwas/pheno/")
)

by_dest <- split(names(out_bucket_map), unlist(out_bucket_map))
for (dest in names(by_dest)) {
  files <- by_dest[[dest]]
  system2("gsutil", c("-m", "cp", files, dest))
}

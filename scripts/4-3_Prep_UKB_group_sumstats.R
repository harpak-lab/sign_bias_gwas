##Preprocessing GWAS sumstats for LDSC (in-house UKB)
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

## input
args <- commandArgs(trailingOnly = TRUE)

if (!("--file" %in% args)) {
  stop("Please specify --file <path>.")
}
if (!("--out" %in% args)) {
  stop("Please specify --out <path>.")
}

file_arg   <- args[which(args=="--file")+1]
out_arg   <- args[which(args=="--out")+1]   

print("Processing...")

## Load gwas
gwas <- fread(file_arg, sep="\t")
gwas <- gwas %>% select("rsid"=ID, "chr"="#CHROM", "pos"=POS, "a1"=A1, "a2"=OMITTED, "beta"=BETA,"se"=SE,"pval"=P)

#export
fwrite(gwas, paste0(out_arg, "_filtered.txt"), sep="\t")

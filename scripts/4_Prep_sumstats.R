##Preprocessing GWAS sumstats for LDSC and S-LD4M
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

## input
args <- commandArgs(trailingOnly = TRUE)

if (!("--file" %in% args)) {
  stop("Please specify --file <path>.")
}
if (!("--snp" %in% args)) {
  stop("Please specify --snp <path>.")
}
if (!("--type" %in% args)) {
  stop("Please specify --type, UKB or AOU.")
}
if (!("--out" %in% args)) {
  stop("Please specify --out <path>.")
}

file_arg   <- args[which(args=="--file")+1]
type_arg   <- args[which(args=="--type")+1]
snp_arg   <- args[which(args=="--snp")+1]   
out_arg   <- args[which(args=="--out")+1]   

print("Processing...")

## Load gwas
gwas <- fread(file_arg, sep="\t")

# load snps
snps <- fread(snp_arg, sep="\t")

# process gwas sumstatf
if (type_arg=="AOU") {
  
  #rename
  names(gwas)[1] <- "CHR"
  gwas[, CHR := paste0("chr", CHR)]

  #as.character
  snps$POS <- as.character(snps$POS)
  gwas$POS <- as.character(gwas$POS)
  
  # merge
  gwas <- merge(gwas, snps, by=c("CHR", "POS"))
  
  # remove duplicates
  gwas_unique_rsid <- gwas[!duplicated(gwas$rsid), ]
  
  gwas_unique_rsid[, A2 := fifelse(A1 == REF, ALT,
                       fifelse(A1 == ALT, REF, NA))]
  
  ## select cols
  gwas_filt <- gwas_unique_rsid %>% select(rsid, chr=CHR, pos=POS, 
                                           Allele1=A1, Allele2=A2, beta=BETA, se=SE, pval=P)
  gwas_filt$chr <- type.convert(sub("^chr", "", gwas_filt$chr), as.is = TRUE)
  gwas_filt$pos <- as.integer(gwas_filt$pos)
  gwas_filt <- gwas_filt[order(gwas_filt$chr, gwas_filt$pos, decreasing = FALSE), ]
  
} else {
  #bind with variants
  gwas <- cbind(gwas,snps)
  
  # remove duplicates
  gwas_unique_rsid <- gwas[!duplicated(gwas$rsid), ]
  
  # keep info > 0.9
  gwas_unique_rsid <- gwas_unique_rsid[info > 0.9]
  
  ## select cols
  gwas_filt <- gwas_unique_rsid %>% select(rsid, chr, pos, 
                                           ref, alt, beta, se, pval)
}

## convert to chisq and drop rs (for use in S-LD4M)
betas <- gwas_filt[,c(1,6,7)]
betas <- betas %>%
  mutate(
    rsid = str_remove(rsid, "^rs")  # remove "rs" prefix; leaves others unchanged
  )

betas$beta <- as.numeric(betas$beta)
betas$se <- as.numeric(betas$se)

betas$chisq <- (betas$beta / betas$se)^2

chisq <- betas[,c(1,4)]

#export
fwrite(gwas_filt, paste0(out_arg, "_filtered.txt"), sep="\t")
fwrite(chisq, paste0(out_arg, "_chisq.txt"), sep="\t")


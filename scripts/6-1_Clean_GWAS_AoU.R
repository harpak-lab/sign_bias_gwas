## Clean AoU GWAS

#libraries
library("dplyr")
library("dtplyr")
library("tidyr")
library("data.table")
library("Biostrings")
library("magrittr")
library("rtracklayer")
library("GenomicRanges")
library("R.utils")

## Read in GWAS
gwas_file <- commandArgs(trailingOnly = TRUE)[1]
out_file  <- commandArgs(trailingOnly = TRUE)[2]

sdat <- fread(gwas_file, sep = "\t")
sdat <- sdat[ERRCODE == "."]

sdat$AC_Allele2 <- (2*sdat$OBS_CT)*sdat$A1_FREQ

sdat <- sdat[,c(1,2,3,4,5,17,9,12,13,15)]

names(sdat) <- c("CHR","POS", "MarkerID", "Allele1","Allele2",
  "AC_Allele2", "AF_Allele2","BETA","SE","Pvalue")

sdat <- sdat[ , if (.N == 1) .SD, by = .(CHR,POS)]
valid <- c("A","C","G","T")

sdat <- sdat[
  nchar(Allele1) == 1 & nchar(Allele2) == 1 &
    Allele1 %in% valid  & Allele2 %in% valid
]

names(sdat)[c(1,2)] <- c("chr","pos")
sdat$chr <- type.convert(sub("^chr", "", sdat$chr), as.is = TRUE)

#edit columns
sdat$pos <- as.numeric(sdat$pos)
sdat$AF_Allele2 <- as.numeric(sdat$AF_Allele2)
sdat$AC_Allele2 <- as.numeric(sdat$AC_Allele2)
sdat$BETA  <- as.numeric(sdat$BETA)
sdat$SE <- as.numeric(sdat$SE)
sdat$Pvalue <- as.numeric(sdat$Pvalue)
aou_rsid <- fread(file="rsids_from_ukb.tsv",
                  sep="\t")
names(aou_rsid) <- c("chr","pos","rsid")
aou_rsid$chr <- sub("^chr", "", aou_rsid$chr)
aou_rsid$chr <- as.numeric(aou_rsid$chr)
sdat <- merge(sdat, aou_rsid, by = c("chr","pos"))
sdat[, `:=`(
  minor_allele = ifelse(AF_Allele2 <= 0.5, Allele2, Allele1),
  minor_AF         = pmin(AF_Allele2, (1-AF_Allele2))   
)]

#Load LD blocks (from MacDonald et al., 2022)
ld <- read.table("LD_EUR_GRCh38.bed", sep="", 
                 header = TRUE, 
                 check.names=FALSE, 
                 stringsAsFactors = FALSE)
ld$chr <- as.numeric(substr(ld$chr, 4, nchar(ld$chr)))
ld$block_id <- 1:nrow(ld)
sdat <- sdat %>% filter(!chr %in% c("X", "Y")) %>%
  mutate(chr = as.numeric(chr),
         pos = as.numeric(pos)) %>%
  filter(as.numeric(chr) %in% 1:22) %>%
  filter(!is.na(SE)) %>% #remove NA SEs
    
  # polarize effects w.r.t to minor alleles (alt is effect allele in UKB)
  mutate(a2_is_minor = AF_Allele2 <= 0.5,
         m_Beta = ifelse(a2_is_minor, BETA, -BETA),
         MAF = minor_AF)
setDT(ld)
sdat <- sdat %>%
  left_join(
    ld,
    join_by(chr, between(pos, start, end))   
  ) %>%
  select(-start, -end)
sdat <- sdat[!is.na(block_id)]

#MAC >= 20
sdat<- sdat[AC_Allele2 >= 20]

#Export cleaned GWAS 
sdat <- sdat %>% select(chr, pos, rsid, MarkerID, Allele1,Allele2,AC_Allele2,AF_Allele2,
                BETA,SE,Pvalue,minor_allele,minor_AF,a2_is_minor,
                m_Beta,MAF,block_id)
fwrite(sdat,
       file      = out_file,   
       sep       = "\t",           # tab delimiter
       quote     = FALSE,          # don’t quote strings
       showProgress = TRUE)

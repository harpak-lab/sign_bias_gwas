library(data.table)
library(dplyr)

#load VAT snps
vat_snps <- fread("vat_rsids.clean.tsv", sep="\t")

#Load Neale UKB variants 
neale_snps <- fread("Neale.variants.tsv", sep="\t")
vat_ukb_subset <- vat_snps[vat_snps$dbsnp_rsid %in% neale_snps$rsid, ]

fwrite(vat_ukb_subset, "rsids_from_ukb.tsv", sep="\t")

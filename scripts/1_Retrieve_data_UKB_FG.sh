#!/bin/bash
############# Get UKB phenotype data and GWAS summmary statistics 
############# (UKB and FinnGen), AoU data queried from Workbench

############# Query UKB phenotype data 
## 31 = reported sex
## 21000 = self-reported ethnicity (instance 0)
## 21022 = age at assessment
## 22000 = batch
## 22001 = genetic sex
## 22003 = heterozygosity
## 22004 = heterozygosity, PCA corrected
## 22005 = missingness
## 22006 = genetic ethnic grouping
## 22019 = sex chromosome aneuploidy
## 22020 = used in genetic PCs
## 22021 = kinship summary
## 22022 = X intensity
## 22023 = Y intensity
## 22024 = DNA concentration
## 22025 = Cluster.CR
## 22026 = dQC
## 22027 = outliers for het/missing
## 22028 = use in phasing chr1-22
## 22029 = use in phasing chrX
## 22030  = use in phasing chrXY

## 50 = standing height
## 21001 = BMI
## 21002 = weight
## 30190 = monocyte percentage
## 30200 = neutrophill percentage
## 30220 = basophil percentage
## 22127 = asthma
## 40001 = underlying (primary) cause of death (IDC10)
## 40002 = contributory (secondary) causes of death
## 41270 = diagnoses - ICD10

cat <<EOF >> pheno_fields_ukb.txt
31 
21000 
21022 
22000 
22001 
22003 
22004 
22005 
22006 
22019 
22020 
22021 
22022 
22023 
22024 
22025 
22026 
22027 
22028 
20029
22030
EOF

cat <<EOF >> qc_fields_ukb.txt
50
21001
21002
30190
30200
30220
22127 
40001
40002
41270
EOF

# Get indicated phenotype values
./ukbconv ukb677255.enc_ukb csv -ipheno_fields_ukb.txt -oUKB.phenos.all.csv
./ukbconv ukb677255.enc_ukb csv -iqc_fields_ukb.txt -oUKB.qc.csv

# Get UKB GWAS variant files (Neale)
wget wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz -O Neale.variants.tsv

# Get UKB GWAS summary statistics (raw units)
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/50_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O Neale.height.tsv.bgz 
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O Neale.BMI.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21002_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O Neale.weight.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30190_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O Neale.monocyte_pct.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30200_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O Neale.neutrophil_pct.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30220_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O Neale.basophil_pct.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/E4_DM1.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.T1D.tsv.gz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/E4_DM2.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.T2D.tsv.gz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/22127.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.asthma.tsv.gz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/F20.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.schizophrenia.tsv.gz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/AD.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.AD.tsv.gz


############ Retrieve FinnGen GWAS summary statistics (r12; google cloud platform)
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_T1D.gz
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_T2D.gz
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_J10_ASTHMA_EXMORE.gz
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_G6_ALZHEIMER.gz
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_F5_SCHZPHR.gz

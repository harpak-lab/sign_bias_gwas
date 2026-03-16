#!/bin/bash

############# Get UKB phenotype data and GWAS summmary statistics 
#(UKB and FinnGen), AoU data queried from Workbench

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
## 22026 = Affymetrix quality control metric "dQC"
## 22027 = outliers for het/missing
## 22028 = use in phasing chr1-22
## 22029 = use in phasing chrX
## 22030  = use in phasing chrXY
## 22009  = genetic PCs

## 50 = standing height
## 21001 = BMI
## 21002 = weight
## 30190 = monocyte percentage
## 30200 = neutrophill percentage
## 30220 = basophil percentage
## 41202 = main ICD10 diagnoses
## 40001 = underlying (primary) cause of death (IDC10)
## 40002 = contributory (secondary) causes of death
## 41270 = diagnoses - ICD10
## 30000 = white blood cell count
## 30010 = red blood cell count
## 30050 = mean corpuscular haemoglobin

cat <<EOF >> qc_fields_ukb.txt
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
22029
22030
22009
EOF

cat <<EOF >> pheno_fields_ukb.txt
50
21001
21002
30190
30200
30220
41202
40001
40002
41270
30000
30010
30050
EOF

# Get indicated phenotype values
./ukbconv ukb677255.enc_ukb csv -ipheno_fields_ukb.txt -oUKB.phenos.all
./ukbconv ukb677255.enc_ukb csv -iqc_fields_ukb.txt -oUKB.qc

# Get UKB GWAS variant files (Neale)
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz -O Neale.variants.tsv

# Get UKB GWAS summary statistics (raw units)
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/50_raw.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.height.tsv.bgz 
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_raw.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.BMI.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21002_raw.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.weight.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30190_raw.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.monocyte_percentage.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30200_raw.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.neutrophil_percentage.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30220_raw.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.basophil_percentage.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30000_raw.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.white_blood_cell_count.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30010_raw.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.red_blood_cell_count.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30050_raw.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.mean_corpuscular_hemoglobin.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/J45.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.asthma_J45.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/E4_DM1.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.type1_diabetes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/E4_DM2.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.type2_diabetes.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/F20.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.schizophrenia_ICD10.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/AD.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.AD.tsv.bgz

# Get UKB GWAS summary statistics (transformed; irnt)
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.height.irnt.tsv.bgz 
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.BMI.irnt.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21002_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.weight.irnt.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30190_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.monocyte_percentage.irnt.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30200_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.neutrophil_percentage.irnt.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30220_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.basophil_percentage.irnt.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30000_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.white_blood_cell_count.irnt.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30010_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.red_blood_cell_count.irnt.tsv.bgz
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30050_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O Neale.mean_corpuscular_hemoglobin.irnt.tsv.bgz

############ Retrieve FinnGen GWAS summary statistics (r12; google cloud platform)
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_T1D.gz
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_T2D.gz
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_J10_ASTHMA_EXMORE.gz
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_G6_ALZHEIMER.gz
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_F5_SCHZPHR.gz

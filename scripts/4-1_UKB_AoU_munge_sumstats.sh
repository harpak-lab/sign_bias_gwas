#!/bin/bash
PYTHON="${PYTHON:-python}"
MUNGE="/software/ldsc/munge_sumstats.py"

IN_DIR="/munged/gwas_processed"
MERGE_ALLELES="/eur_w_ld_chr/w_hm3.snplist"

CHUNKSIZE="${CHUNKSIZE:-500000}"

A1="Allele2"
A2="Allele1"
SNP="rsid"
P="pval"

run_munge() {
  local sumstats="$1"
  local outprefix="$2"
  shift 2
  echo "==> munging: $(basename "$sumstats") -> ${outprefix}.sumstats.gz"
  "$PYTHON" "$MUNGE" \
    --sumstats "$sumstats" \
    --a1 "$A1" --a2 "$A2" \
    --snp "$SNP" --p "$P" \
    --out "$outprefix" \
    --merge-alleles "$MERGE_ALLELES" \
    --chunksize "$CHUNKSIZE" \
    "$@"
}

# AoU
run_munge "${IN_DIR}/Asthma_AoU_filtered.txt"               "Asthma_AoU.munge"               --N-cas 9838  --N-con 83119
run_munge "${IN_DIR}/Basophil_percentage_AoU_filtered.txt"  "Basophil_percentage_AoU.munge"  --N 53123
run_munge "${IN_DIR}/BMI_AoU_filtered.txt"                  "BMI_AoU.munge"                  --N 111482
run_munge "${IN_DIR}/Height_AoU_filtered.txt"               "Height_AoU.munge"               --N 111755
run_munge "${IN_DIR}/Monocyte_percentage_AoU_filtered.txt"  "Monocyte_percentage_AoU.munge"  --N 53661
run_munge "${IN_DIR}/Neutrophil_percentage_AoU_filtered.txt""Neutrophil_percentage_AoU.munge" --N 45276
run_munge "${IN_DIR}/Schizophrenia_AoU_filtered.txt"        "Schizophrenia_AoU.munge"        --N-cas 472   --N-con 98864
run_munge "${IN_DIR}/Type_1_diabetes_AoU_filtered.txt"      "Type_1_diabetes_AoU.munge"      --N-cas 1153  --N-con 97447
run_munge "${IN_DIR}/Type_2_diabetes_AoU_filtered.txt"      "Type_2_diabetes_AoU.munge"      --N-cas 10931 --N-con 83621
run_munge "${IN_DIR}/Weight_AoU_filtered.txt"               "Weight_AoU.munge"               --N 111568


# UKB / Neale 
UKB_A1="alt"
UKB_A2="ref"

run_munge "${IN_DIR}/Neale.asthma.txt"                            "Neale.asthma.munge"  "$UKB_A1" "$UKB_A2" --N-cas 11717 --N-con 80070
run_munge "${IN_DIR}/Neale.basophil_percentage_filtered.txt"      "Neale.basophil_percentage.munge"      "$UKB_A1" "$UKB_A2" --N 349861
run_munge "${IN_DIR}/Neale.BMI_filtered.txt"                      "Neale.BMI.munge"                      "$UKB_A1" "$UKB_A2" --N 359983
run_munge "${IN_DIR}/Neale.monocyte_percentage_filtered.txt"      "Neale.monocyte_percentage.munge"      "$UKB_A1" "$UKB_A2" --N 349861
run_munge "${IN_DIR}/Neale.neutrophill_percentage_filtered.txt"   "Neale.neutrophill_percentage.munge"   "$UKB_A1" "$UKB_A2" --N 349861
run_munge "${IN_DIR}/Neale.schizophrenia_ICD10_filtered.txt"      "Neale.schizophrenia_ICD10.munge"      "$UKB_A1" "$UKB_A2" --N-cas 198 --N-con 360996
run_munge "${IN_DIR}/Neale.standing_height_filtered.txt"          "Neale.standing_height.munge"          "$UKB_A1" "$UKB_A2" --N 360388
run_munge "${IN_DIR}/Neale.type1_diabetes_filtered.txt"           "Neale.type1_diabetes.munge"           "$UKB_A1" "$UKB_A2" --N-cas 583 --N-con 360611
run_munge "${IN_DIR}/Neale.type2_diabetes_filtered.txt"           "Neale.type2_diabetes.munge"           "$UKB_A1" "$UKB_A2" --N-cas 888 --N-con 360306
run_munge "${IN_DIR}/Neale.weight_filtered.txt"                   "Neale.weight.munge"                   "$UKB_A1" "$UKB_A2" --N 354838


echo "All munge jobs completed."

#!/bin/bash

#input
PYTHON="${PYTHON:-python}"
RSCRIPT="${RSCRIPT:-Rscript}"

PREP_SCRIPT="/scripts/4_Prep_sumstats.R"
MUNGE="/software/ldsc/munge_sumstats.py"

IN_DIR="/munged/gwas_raw"
OUT_DIR="/munged/gwas_processed"
MERGE_ALLELES="/eur_w_ld_chr/w_hm3.snplist"

CHUNKSIZE="${CHUNKSIZE:-500000}"

mkdir -p "$OUT_DIR"

#run functions
run_prep() {
  local infile="$1"
  local snpfile="$2"
  local type="$3"
  local outprefix="$4"

  echo "preprocessing: $(basename "$infile") -> ${outprefix}_filtered.txt"
  "$RSCRIPT" "$PREP_SCRIPT" \
    --file "$infile" \
    --snp "$snpfile" \
    --type "$type" \
    --out "$outprefix"
}

run_munge() {
  local sumstats="$1"
  local outprefix="$2"
  local a1_col="$3"
  local a2_col="$4"
  shift 4

  echo "munging: $(basename "$sumstats") -> ${outprefix}.sumstats.gz"
  "$PYTHON" "$MUNGE" \
    --sumstats "$sumstats" \
    --a1 "$a1_col" \
    --a2 "$a2_col" \
    --snp rsid \
    --p pval \
    --out "$outprefix" \
    --merge-alleles "$MERGE_ALLELES" \
    --chunksize "$CHUNKSIZE" \
    "$@"
}

# SNPs
AOU_SNP="/data/Neale.variants.tsv"
UKB_SNP="/hm3/hm3_snplist_hg38.txt"

# AoU: preprocess then munge
AOU_A1="Allele1"
AOU_A2="Allele2"

#run prep and munge (AOU)
run_prep "${IN_DIR}/Alzheimers_AoU.txt"              "$AOU_SNP" "AOU" "${OUT_DIR}/Alzheimers_AoU"
run_munge "${OUT_DIR}/Alzheimers_AoU_filtered.txt"   "${OUT_DIR}/Alzheimers_AoU.munge" "$AOU_A1" "$AOU_A2" --N-cas 367 --N-con 161007

run_prep "${IN_DIR}/Asthma_AoU.txt"                  "$AOU_SNP" "AOU" "${OUT_DIR}/Asthma_AoU"
run_munge "${OUT_DIR}/Asthma_AoU_filtered.txt"       "${OUT_DIR}/Asthma_AoU.munge" "$AOU_A1" "$AOU_A2" --N-cas 19910 --N-con 134877

run_prep "${IN_DIR}/Basophil_percentage_AoU.txt"     "$AOU_SNP" "AOU" "${OUT_DIR}/Basophil_percentage_AoU"
run_munge "${OUT_DIR}/Basophil_percentage_AoU_filtered.txt" \
          "${OUT_DIR}/Basophil_percentage_AoU.munge" "$AOU_A1" "$AOU_A2" --N 78988

run_prep "${IN_DIR}/BMI_AoU.txt"                     "$AOU_SNP" "AOU" "${OUT_DIR}/BMI_AoU"
run_munge "${OUT_DIR}/BMI_AoU_filtered.txt"          "${OUT_DIR}/BMI_AoU.munge" "$AOU_A1" "$AOU_A2" --N 198762

run_prep "${IN_DIR}/Height_AoU.txt"                  "$AOU_SNP" "AOU" "${OUT_DIR}/Height_AoU"
run_munge "${OUT_DIR}/Height_AoU_filtered.txt"       "${OUT_DIR}/Height_AoU.munge" "$AOU_A1" "$AOU_A2" --N 199529

run_prep "${IN_DIR}/Monocyte_percentage_AoU.txt"     "$AOU_SNP" "AOU" "${OUT_DIR}/Monocyte_percentage_AoU"
run_munge "${OUT_DIR}/Monocyte_percentage_AoU_filtered.txt" \
          "${OUT_DIR}/Monocyte_percentage_AoU.munge" "$AOU_A1" "$AOU_A2" --N 79163

run_prep "${IN_DIR}/Neutrophil_percentage_AoU.txt"   "$AOU_SNP" "AOU" "${OUT_DIR}/Neutrophil_percentage_AoU"
run_munge "${OUT_DIR}/Neutrophil_percentage_AoU_filtered.txt" \
          "${OUT_DIR}/Neutrophil_percentage_AoU.munge" "$AOU_A1" "$AOU_A2" --N 69249

run_prep "${IN_DIR}/Schizophrenia_AoU.txt"           "$AOU_SNP" "AOU" "${OUT_DIR}/Schizophrenia_AoU"
run_munge "${OUT_DIR}/Schizophrenia_AoU_filtered.txt" \
          "${OUT_DIR}/Schizophrenia_AoU.munge" "$AOU_A1" "$AOU_A2" --N-cas 848 --N-con 160271

run_prep "${IN_DIR}/Type_1_diabetes_AoU.txt"         "$AOU_SNP" "AOU" "${OUT_DIR}/Type_1_diabetes_AoU"
run_munge "${OUT_DIR}/Type_1_diabetes_AoU_filtered.txt" \
          "${OUT_DIR}/Type_1_diabetes_AoU.munge" "$AOU_A1" "$AOU_A2" --N-cas 2263 --N-con 158166

run_prep "${IN_DIR}/Type_2_diabetes_AoU.txt"         "$AOU_SNP" "AOU" "${OUT_DIR}/Type_2_diabetes_AoU"
run_munge "${OUT_DIR}/Type_2_diabetes_AoU_filtered.txt" \
          "${OUT_DIR}/Type_2_diabetes_AoU.munge" "$AOU_A1" "$AOU_A2" --N-cas 22805 --N-con 134033

run_prep "${IN_DIR}/Weight_AoU.txt"                  "$AOU_SNP" "AOU" "${OUT_DIR}/Weight_AoU"
run_munge "${OUT_DIR}/Weight_AoU_filtered.txt"       "${OUT_DIR}/Weight_AoU.munge" "$AOU_A1" "$AOU_A2" --N 199194

run_prep "${IN_DIR}/White_blood_cell_count_AoU.txt"  "$AOU_SNP" "AOU" "${OUT_DIR}/White_blood_cell_count_AoU"
run_munge "${OUT_DIR}/White_blood_cell_count_AoU_filtered.txt" \
          "${OUT_DIR}/White_blood_cell_count_AoU.munge" "$AOU_A1" "$AOU_A2" --N 71892

run_prep "${IN_DIR}/Red_blood_cell_count_AoU.txt"    "$AOU_SNP" "AOU" "${OUT_DIR}/Red_blood_cell_count_AoU"
run_munge "${OUT_DIR}/Red_blood_cell_count_AoU_filtered.txt" \
          "${OUT_DIR}/Red_blood_cell_count_AoU.munge" "$AOU_A1" "$AOU_A2" --N 68492

run_prep "${IN_DIR}/Mean_corpuscular_hemoglobin_AoU.txt" "$AOU_SNP" "AOU" "${OUT_DIR}/Mean_corpuscular_hemoglobin_AoU"
run_munge "${OUT_DIR}/Mean_corpuscular_hemoglobin_AoU_filtered.txt" \
          "${OUT_DIR}/Mean_corpuscular_hemoglobin_AoU.munge" "$AOU_A1" "$AOU_A2" --N 96648


# UKB / Neale, prep and munge

UKB_A1="alt"
UKB_A2="ref"

run_prep "${IN_DIR}/Neale.alzheimers.txt"                    "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.alzheimers"
run_munge "${OUT_DIR}/Neale.alzheimers_filtered.txt"        "${OUT_DIR}/Neale.alzheimers.munge" "$UKB_A1" "$UKB_A2" --N-cas 119 --N-con 361075

run_prep "${IN_DIR}/Neale.asthma.txt"                       "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.asthma"
run_munge "${OUT_DIR}/Neale.asthma_filtered.txt"           "${OUT_DIR}/Neale.asthma.munge" "$UKB_A1" "$UKB_A2" --N-cas 1693 --N-con 359501

run_prep "${IN_DIR}/Neale.basophil_percentage.txt"          "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.basophil_percentage"
run_munge "${OUT_DIR}/Neale.basophil_percentage_filtered.txt" \
          "${OUT_DIR}/Neale.basophil_percentage.munge" "$UKB_A1" "$UKB_A2" --N 349861

run_prep "${IN_DIR}/Neale.BMI.txt"                          "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.BMI"
run_munge "${OUT_DIR}/Neale.BMI_filtered.txt"              "${OUT_DIR}/Neale.BMI.munge" "$UKB_A1" "$UKB_A2" --N 359983

run_prep "${IN_DIR}/Neale.monocyte_percentage.txt"         "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.monocyte_percentage"
run_munge "${OUT_DIR}/Neale.monocyte_percentage_filtered.txt" \
          "${OUT_DIR}/Neale.monocyte_percentage.munge" "$UKB_A1" "$UKB_A2" --N 349861

run_prep "${IN_DIR}/Neale.neutrophil_percentage.txt"       "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.neutrophil_percentage"
run_munge "${OUT_DIR}/Neale.neutrophil_percentage_filtered.txt" \
          "${OUT_DIR}/Neale.neutrophil_percentage.munge" "$UKB_A1" "$UKB_A2" --N 349861

run_prep "${IN_DIR}/Neale.schizophrenia_ICD10.txt"         "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.schizophrenia_ICD10"
run_munge "${OUT_DIR}/Neale.schizophrenia_ICD10_filtered.txt" \
          "${OUT_DIR}/Neale.schizophrenia_ICD10.munge" "$UKB_A1" "$UKB_A2" --N-cas 198 --N-con 360996

run_prep "${IN_DIR}/Neale.standing_height.txt"             "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.standing_height"
run_munge "${OUT_DIR}/Neale.standing_height_filtered.txt" \
          "${OUT_DIR}/Neale.standing_height.munge" "$UKB_A1" "$UKB_A2" --N 360388

run_prep "${IN_DIR}/Neale.type1_diabetes.txt"              "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.type1_diabetes"
run_munge "${OUT_DIR}/Neale.type1_diabetes_filtered.txt" \
          "${OUT_DIR}/Neale.type1_diabetes.munge" "$UKB_A1" "$UKB_A2" --N-cas 583 --N-con 360611

run_prep "${IN_DIR}/Neale.type2_diabetes.txt"              "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.type2_diabetes"
run_munge "${OUT_DIR}/Neale.type2_diabetes_filtered.txt" \
          "${OUT_DIR}/Neale.type2_diabetes.munge" "$UKB_A1" "$UKB_A2" --N-cas 888 --N-con 360306

run_prep "${IN_DIR}/Neale.weight.txt"                      "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.weight"
run_munge "${OUT_DIR}/Neale.weight_filtered.txt"          "${OUT_DIR}/Neale.weight.munge" "$UKB_A1" "$UKB_A2" --N 354838

run_prep "${IN_DIR}/Neale.white_blood_cell_count.txt"      "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.white_blood_cell_count"
run_munge "${OUT_DIR}/Neale.white_blood_cell_count_filtered.txt" \
          "${OUT_DIR}/Neale.white_blood_cell_count.munge" "$UKB_A1" "$UKB_A2" --N 350470

run_prep "${IN_DIR}/Neale.red_blood_cell_count.txt"        "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.red_blood_cell_count"
run_munge "${OUT_DIR}/Neale.red_blood_cell_count_filtered.txt" \
          "${OUT_DIR}/Neale.red_blood_cell_count.munge" "$UKB_A1" "$UKB_A2" --N 350475

run_prep "${IN_DIR}/Neale.mean_corpuscular_hemoglobin.txt" "$UKB_SNP" "UKB" "${OUT_DIR}/Neale.mean_corpuscular_hemoglobin"
run_munge "${OUT_DIR}/Neale.mean_corpuscular_hemoglobin_filtered.txt" \
          "${OUT_DIR}/Neale.mean_corpuscular_hemoglobin.munge" "$UKB_A1" "$UKB_A2" --N 350472
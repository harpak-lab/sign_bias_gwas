#!/bin/bash

### Run cleaning steps

RSCRIPT="${RSCRIPT:-Rscript}"

# Scripts
AOU_SCRIPT="/stor/work/Kirkpatrick/scratch/Jared/signbias/UKB_round2/scripts/6-1_Clean_GWAS_AoU.R"
UKB_SCRIPT="/stor/work/Kirkpatrick/scratch/Jared/signbias/UKB_round2/scripts/6-2_Clean_GWAS_UKB.R"
FG_SCRIPT="/stor/work/Kirkpatrick/scratch/Jared/signbias/UKB_round2/scripts/6-3_Clean_GWAS_FG.R"

# LD blocks bed (used by UKB script)
LD_BED="/stor/work/Kirkpatrick/scratch/Jared/signbias/data/LD_EUR.bed"

# Output directories
OUT_AOU="/stor/work/Kirkpatrick/scratch/Jared/signbias/AoU_round2/cleaned3"
OUT_UKB="/stor/work/Kirkpatrick/scratch/Jared/signbias/UKB_round2/cleaned3"
OUT_FG="/stor/work/Kirkpatrick/scratch/Jared/signbias/FG_round2/cleaned3"

LOG_AOU="${OUT_AOU}/logs"
LOG_UKB="${OUT_UKB}/logs"
LOG_FG="${OUT_FG}/logs"

mkdir -p "$OUT_AOU" "$OUT_UKB" "$OUT_FG" "$LOG_AOU" "$LOG_UKB" "$LOG_FG"

run_aou () {
  local in_file="$1"
  local out_file="$2"
  local log_file="$3"

  echo "==> AoU clean: $(basename "$in_file")"
  "$RSCRIPT" "$AOU_SCRIPT" "$in_file" "$out_file" >> "$log_file" 2>&1
}

run_ukb () {
  local in_file="$1"
  local out_file="$2"
  local log_file="$3"

  echo "==> UKB clean: $(basename "$in_file")"
  "$RSCRIPT" "$UKB_SCRIPT" -f "$in_file" -b "$LD_BED" -o "$out_file" >> "$log_file" 2>&1
}

run_fg () {
  local in_file="$1"
  local out_file="$2"
  local log_file="$3"

  echo "==> FinnGen clean: $(basename "$in_file")"
  "$RSCRIPT" "$FG_SCRIPT" "$in_file" "$out_file" >> "$log_file" 2>&1
}

UKB_IN_DIR="/stor/work/Kirkpatrick/scratch/Jared/signbias/data/GWAS/Neale/round2"

run_ukb "${UKB_IN_DIR}/Neale.basophil_percentage.tsv.bgz"    "${OUT_UKB}/Neale.basophil_percentage.cleaned.txt"     "${LOG_UKB}/Neale.basophil_percentage.log"
run_ukb "${UKB_IN_DIR}/Neale.BMI.tsv.bgz"                    "${OUT_UKB}/Neale.BMI.cleaned.txt"                     "${LOG_UKB}/Neale.BMI.log"
run_ukb "${UKB_IN_DIR}/Neale.monocyte_percentage.tsv.bgz"    "${OUT_UKB}/Neale.monocyte_percentage.cleaned.txt"     "${LOG_UKB}/Neale.monocyte_percentage.log"
run_ukb "${UKB_IN_DIR}/Neale.neutrophill_percentage.tsv.bgz" "${OUT_UKB}/Neale.neutrophill_percentage.cleaned.txt"  "${LOG_UKB}/Neale.neutrophill_percentage.log"
run_ukb "${UKB_IN_DIR}/Neale.schizophrenia_ICD10.tsv.bgz"    "${OUT_UKB}/Neale.schizophrenia_ICD10.cleaned.txt"     "${LOG_UKB}/Neale.schizophrenia_ICD10.log"
run_ukb "${UKB_IN_DIR}/Neale.standing_height.tsv.bgz"        "${OUT_UKB}/Neale.standing_height.cleaned.txt"         "${LOG_UKB}/Neale.standing_height.log"
run_ukb "${UKB_IN_DIR}/Neale.type1_diabetes.tsv.bgz"         "${OUT_UKB}/Neale.type1_diabetes.cleaned.txt"          "${LOG_UKB}/Neale.type1_diabetes.log"
run_ukb "${UKB_IN_DIR}/Neale.type2_diabetes.tsv.bgz"         "${OUT_UKB}/Neale.type2_diabetes.cleaned.txt"          "${LOG_UKB}/Neale.type2_diabetes.log"
run_ukb "${UKB_IN_DIR}/Neale.weight.tsv.bgz"                 "${OUT_UKB}/Neale.weight.cleaned.txt"                  "${LOG_UKB}/Neale.weight.log"

AOU_IN_DIR="/stor/work/Kirkpatrick/scratch/Jared/signbias/data/GWAS/AoU/round2"
run_aou "${AOU_IN_DIR}/Asthma_AoU.tsv.gz"                "${OUT_AOU}/Asthma_AoU.cleaned.txt"                "${LOG_AOU}/Asthma_AoU.log"
run_aou "${AOU_IN_DIR}/Basophil_percentage_AoU.tsv.gz"   "${OUT_AOU}/Basophil_percentage_AoU.cleaned.txt"   "${LOG_AOU}/Basophil_percentage_AoU.log"
run_aou "${AOU_IN_DIR}/BMI_AoU.tsv.gz"                   "${OUT_AOU}/BMI_AoU.cleaned.txt"                   "${LOG_AOU}/BMI_AoU.log"
run_aou "${AOU_IN_DIR}/Height_AoU.tsv.gz"                "${OUT_AOU}/Height_AoU.cleaned.txt"                "${LOG_AOU}/Height_AoU.log"
run_aou "${AOU_IN_DIR}/Monocyte_percentage_AoU.tsv.gz"   "${OUT_AOU}/Monocyte_percentage_AoU.cleaned.txt"   "${LOG_AOU}/Monocyte_percentage_AoU.log"
run_aou "${AOU_IN_DIR}/Neutrophil_percentage_AoU.tsv.gz" "${OUT_AOU}/Neutrophil_percentage_AoU.cleaned.txt" "${LOG_AOU}/Neutrophil_percentage_AoU.log"
run_aou "${AOU_IN_DIR}/Schizophrenia_AoU.tsv.gz"         "${OUT_AOU}/Schizophrenia_AoU.cleaned.txt"         "${LOG_AOU}/Schizophrenia_AoU.log"
run_aou "${AOU_IN_DIR}/Type_1_diabetes_AoU.tsv.gz"       "${OUT_AOU}/Type_1_diabetes_AoU.cleaned.txt"       "${LOG_AOU}/Type_1_diabetes_AoU.log"
run_aou "${AOU_IN_DIR}/Type_2_diabetes_AoU.tsv.gz"       "${OUT_AOU}/Type_2_diabetes_AoU.cleaned.txt"       "${LOG_AOU}/Type_2_diabetes_AoU.log"
run_aou "${AOU_IN_DIR}/Weight_AoU.tsv.gz"                "${OUT_AOU}/Weight_AoU.cleaned.txt"                "${LOG_AOU}/Weight_AoU.log"

FG_IN_DIR="/stor/work/Kirkpatrick/scratch/Jared/signbias/data/GWAS/FinnGen/R12"
run_fg "${FG_IN_DIR}/finngen_R12_F5_SCHZPHR.tsv.gz"          "${OUT_FG}/finngen_R12_F5_SCHZPHR.cleaned.txt"          "${LOG_FG}/finngen_R12_F5_SCHZPHR.log"
run_fg "${FG_IN_DIR}/finngen_R12_G6_ALZHEIMER.tsv.gz"        "${OUT_FG}/finngen_R12_G6_ALZHEIMER.cleaned.txt"        "${LOG_FG}/finngen_R12_G6_ALZHEIMER.log"
run_fg "${FG_IN_DIR}/finngen_R12_J10_ASTHMA_EXMORE.tsv.gz"   "${OUT_FG}/finngen_R12_J10_ASTHMA_EXMORE.cleaned.txt"   "${LOG_FG}/finngen_R12_J10_ASTHMA_EXMORE.log"
run_fg "${FG_IN_DIR}/finngen_R12_T1D.tsv.gz"                 "${OUT_FG}/finngen_R12_T1D.cleaned.txt"                 "${LOG_FG}/finngen_R12_T1D.log"
run_fg "${FG_IN_DIR}/finngen_R12_T2D.tsv.gz"                 "${OUT_FG}/finngen_R12_T2D.cleaned.txt"                 "${LOG_FG}/finngen_R12_T2D.log"

echo "Done."
echo "UKB outputs: $OUT_UKB"
echo "AoU outputs: $OUT_AOU"
echo "FinnGen outputs: $OUT_FG"



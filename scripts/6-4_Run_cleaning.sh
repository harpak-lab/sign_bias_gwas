#!/bin/bash

### Run cleaning steps

RSCRIPT="${RSCRIPT:-Rscript}"

# Scripts
AOU_SCRIPT="/scripts/6-1_Clean_GWAS_AoU.R"
UKB_SCRIPT="/scripts/6-2_Clean_GWAS_UKB.R"
FG_SCRIPT="/scripts/6-3_Clean_GWAS_FG.R"

# LD blocks bed (used by UKB script)
LD_BED="/data/LD_EUR.bed"

# Output directories
OUT_AOU="/cleaned"
OUT_UKB="/cleaned"
OUT_FG="/cleaned"

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

UKB_IN_DIR="/GWAS/Neale/"
run_ukb "${UKB_IN_DIR}/Neale.AD.tsv.bgz"                          "${OUT_UKB}/Neale.alzheimers.cleaned.txt"                   "${LOG_UKB}/Neale.alzheimers.log"
run_ukb "${UKB_IN_DIR}/Neale.asthma_J45.tsv.bgz"                  "${OUT_UKB}/Neale.asthma.cleaned.txt"                       "${LOG_UKB}/Neale.asthma.log"
run_ukb "${UKB_IN_DIR}/Neale.basophil_percentage.tsv.bgz"         "${OUT_UKB}/Neale.basophil_percentage.cleaned.txt"          "${LOG_UKB}/Neale.basophil_percentage.log"
run_ukb "${UKB_IN_DIR}/Neale.BMI.tsv.bgz"                         "${OUT_UKB}/Neale.BMI.cleaned.txt"                          "${LOG_UKB}/Neale.BMI.log"
run_ukb "${UKB_IN_DIR}/Neale.monocyte_percentage.tsv.bgz"         "${OUT_UKB}/Neale.monocyte_percentage.cleaned.txt"          "${LOG_UKB}/Neale.monocyte_percentage.log"
run_ukb "${UKB_IN_DIR}/Neale.neutrophill_percentage.tsv.bgz"      "${OUT_UKB}/Neale.neutrophill_percentage.cleaned.txt"       "${LOG_UKB}/Neale.neutrophill_percentage.log"
run_ukb "${UKB_IN_DIR}/Neale.white_blood_cell_count.tsv.bgz"      "${OUT_UKB}/Neale.white_blood_cell_count.cleaned.txt"       "${LOG_UKB}/Neale.white_blood_cell_count.log"
run_ukb "${UKB_IN_DIR}/Neale.red_blood_cell_count.tsv.bgz"        "${OUT_UKB}/Neale.red_blood_cell_count.cleaned.txt"         "${LOG_UKB}/Neale.red_blood_cell_count.log"
run_ukb "${UKB_IN_DIR}/Neale.mean_corpuscular_hemoglobin.tsv.bgz" "${OUT_UKB}/Neale.mean_corpuscular_hemoglobin.cleaned.txt"  "${LOG_UKB}/Neale.mean_corpuscular_hemoglobin.log"
run_ukb "${UKB_IN_DIR}/Neale.schizophrenia_ICD10.tsv.bgz"         "${OUT_UKB}/Neale.schizophrenia_ICD10.cleaned.txt"          "${LOG_UKB}/Neale.schizophrenia_ICD10.log"
run_ukb "${UKB_IN_DIR}/Neale.height.tsv.bgz"                      "${OUT_UKB}/Neale.height.cleaned.txt"                       "${LOG_UKB}/Neale.height.log"
run_ukb "${UKB_IN_DIR}/Neale.type1_diabetes.tsv.bgz"              "${OUT_UKB}/Neale.type1_diabetes.cleaned.txt"               "${LOG_UKB}/Neale.type1_diabetes.log"
run_ukb "${UKB_IN_DIR}/Neale.type2_diabetes.tsv.bgz"              "${OUT_UKB}/Neale.type2_diabetes.cleaned.txt"               "${LOG_UKB}/Neale.type2_diabetes.log"
run_ukb "${UKB_IN_DIR}/Neale.weight.tsv.bgz"                      "${OUT_UKB}/Neale.weight.cleaned.txt"                       "${LOG_UKB}/Neale.weight.log"

run_ukb "${UKB_IN_DIR}/Neale.height.irnt.tsv.bgz"                      "${OUT_UKB}/Neale.height.irnt.cleaned.txt"                            "${LOG_UKB}/Neale.height.irnt.log"
run_ukb "${UKB_IN_DIR}/Neale.weight.irnt.tsv.bgz"                      "${OUT_UKB}/Neale.weight.irnt.cleaned.txt"                            "${LOG_UKB}/Neale.weight.irnt.log"
run_ukb "${UKB_IN_DIR}/Neale.BMI.irnt.tsv.bgz"                         "${OUT_UKB}/Neale.BMI.irnt.cleaned.txt"                               "${LOG_UKB}/Neale.BMI.irnt.log"
run_ukb "${UKB_IN_DIR}/Neale.monocyte_percentage.irnt.tsv.bgz"         "${OUT_UKB}/Neale.monocyte_percentage.irnt.cleaned.txt"               "${LOG_UKB}/Neale.monocyte_percentage.irnt.log"
run_ukb "${UKB_IN_DIR}/Neale.basophil_percentage.irnt.tsv.bgz"         "${OUT_UKB}/Neale.basophil_percentage.irnt.cleaned.txt"               "${LOG_UKB}/Neale.basophil_percentage.irnt.log"
run_ukb "${UKB_IN_DIR}/Neale.neutrophill_percentage.irnt.tsv.bgz"      "${OUT_UKB}/Neale.neutrophill_percentage.irnt.cleaned.txt"            "${LOG_UKB}/Neale.neutrophill_percentage.irnt.log"
run_ukb "${UKB_IN_DIR}/Neale.white_blood_cell_count.irnt.tsv.bgz"      "${OUT_UKB}/Neale.white_blood_cell_count.irnt.cleaned.txt"            "${LOG_UKB}/Neale.white_blood_cell_count.irnt.log"
run_ukb "${UKB_IN_DIR}/Neale.red_blood_cell_count.irnt.tsv.bgz"        "${OUT_UKB}/Neale.red_blood_cell_count.irnt.cleaned.txt"              "${LOG_UKB}/Neale.red_blood_cell_count.irnt.log"
run_ukb "${UKB_IN_DIR}/Neale.mean_corpuscular_hemoglobin.irnt.tsv.bgz" "${OUT_UKB}/Neale.mean_corpuscular_hemoglobin.irnt.cleaned.txt"       "${LOG_UKB}/Neale.mean_corpuscular_hemoglobin.irnt.log"

AOU_IN_DIR="/GWAS/AoU/"
run_aou "${AOU_IN_DIR}/Alzheimers_AoU.tsv.gz"                    "${OUT_AOU}/Alzheimers_AoU.cleaned.txt"                     "${LOG_AOU}/Alzheimers_AoU.log"
run_aou "${AOU_IN_DIR}/Asthma_AoU.tsv.gz"                        "${OUT_AOU}/Asthma_AoU.cleaned.txt"                        "${LOG_AOU}/Asthma_AoU.log"
run_aou "${AOU_IN_DIR}/Basophil_percentage_AoU.tsv.gz"           "${OUT_AOU}/Basophil_percentage_AoU.cleaned.txt"           "${LOG_AOU}/Basophil_percentage_AoU.log"
run_aou "${AOU_IN_DIR}/BMI_AoU.tsv.gz"                           "${OUT_AOU}/BMI_AoU.cleaned.txt"                           "${LOG_AOU}/BMI_AoU.log"
run_aou "${AOU_IN_DIR}/Height_AoU.tsv.gz"                        "${OUT_AOU}/Height_AoU.cleaned.txt"                        "${LOG_AOU}/Height_AoU.log"
run_aou "${AOU_IN_DIR}/Monocyte_percentage_AoU.tsv.gz"           "${OUT_AOU}/Monocyte_percentage_AoU.cleaned.txt"           "${LOG_AOU}/Monocyte_percentage_AoU.log"
run_aou "${AOU_IN_DIR}/Neutrophil_percentage_AoU.tsv.gz"         "${OUT_AOU}/Neutrophil_percentage_AoU.cleaned.txt"         "${LOG_AOU}/Neutrophil_percentage_AoU.log"
run_aou "${AOU_IN_DIR}/White_blood_cell_count_AoU.tsv.gz"        "${OUT_AOU}/White_blood_cell_count_AoU.cleaned.txt"        "${LOG_AOU}/White_blood_cell_count_AoU.log"
run_aou "${AOU_IN_DIR}/Red_blood_cell_count_AoU.tsv.gz"          "${OUT_AOU}/Red_blood_cell_count_AoU.cleaned.txt"          "${LOG_AOU}/Red_blood_cell_count_AoU.log"
run_aou "${AOU_IN_DIR}/Mean_corpuscular_hemoglobin_AoU.tsv.gz"   "${OUT_AOU}/Mean_corpuscular_hemoglobin_AoU.cleaned.txt"   "${LOG_AOU}/Mean_corpuscular_hemoglobin_AoU.log"
run_aou "${AOU_IN_DIR}/Schizophrenia_AoU.tsv.gz"                 "${OUT_AOU}/Schizophrenia_AoU.cleaned.txt"                 "${LOG_AOU}/Schizophrenia_AoU.log"
run_aou "${AOU_IN_DIR}/Type_1_diabetes_AoU.tsv.gz"               "${OUT_AOU}/Type_1_diabetes_AoU.cleaned.txt"               "${LOG_AOU}/Type_1_diabetes_AoU.log"
run_aou "${AOU_IN_DIR}/Type_2_diabetes_AoU.tsv.gz"               "${OUT_AOU}/Type_2_diabetes_AoU.cleaned.txt"               "${LOG_AOU}/Type_2_diabetes_AoU.log"
run_aou "${AOU_IN_DIR}/Weight_AoU.tsv.gz"                        "${OUT_AOU}/Weight_AoU.cleaned.txt"                        "${LOG_AOU}/Weight_AoU.log"

FG_IN_DIR="/GWAS/FG"
run_fg "${FG_IN_DIR}/finngen_R12_F5_SCHZPHR.tsv.gz"          "${OUT_FG}/finngen_R12_F5_SCHZPHR.cleaned.txt"          "${LOG_FG}/finngen_R12_F5_SCHZPHR.log"
run_fg "${FG_IN_DIR}/finngen_R12_G6_ALZHEIMER.tsv.gz"        "${OUT_FG}/finngen_R12_G6_ALZHEIMER.cleaned.txt"        "${LOG_FG}/finngen_R12_G6_ALZHEIMER.log"
run_fg "${FG_IN_DIR}/finngen_R12_J10_ASTHMA_EXMORE.tsv.gz"   "${OUT_FG}/finngen_R12_J10_ASTHMA_EXMORE.cleaned.txt"   "${LOG_FG}/finngen_R12_J10_ASTHMA_EXMORE.log"
run_fg "${FG_IN_DIR}/finngen_R12_T1D.tsv.gz"                 "${OUT_FG}/finngen_R12_T1D.cleaned.txt"                 "${LOG_FG}/finngen_R12_T1D.log"
run_fg "${FG_IN_DIR}/finngen_R12_T2D.tsv.gz"                 "${OUT_FG}/finngen_R12_T2D.cleaned.txt"                 "${LOG_FG}/finngen_R12_T2D.log"

echo "Done."
echo "UKB outputs: $OUT_UKB"
echo "AoU outputs: $OUT_AOU"
echo "FinnGen outputs: $OUT_FG"



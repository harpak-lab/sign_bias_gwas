#!/bin/bash

PYTHON="${PYTHON:-python}"
LDSC="/software/ldsc/ldsc.py"

MUNGED_DIR="/munged/"
LDCHR_DIR="/eur_w_ld_chr/"

# where outputs go
OUT_DIR="${OUT_DIR:-ldsc_rg_out}"
mkdir -p "$OUT_DIR"

run_rg () {
  local trait="$1"
  local aou="$2"
  local ukb="$3"
  shift 3

  echo "LDSC: ${trait}"
  "$PYTHON" "$LDSC" \
    --rg "${aou},${ukb}" \
    --ref-ld-chr "$LDCHR_DIR" \
    --w-ld-chr "$LDCHR_DIR" \
    --out "${OUT_DIR}/${trait}" \
    "$@"
}


# Quantitative traits
run_rg "Basophil_percentage" \
  "${MUNGED_DIR}/Basophil_percentage_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.basophil_percentage.munge.sumstats.gz"

run_rg "BMI" \
  "${MUNGED_DIR}/BMI_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.BMI.munge.sumstats.gz"

run_rg "Height" \
  "${MUNGED_DIR}/Height_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.standing_height.munge.sumstats.gz"

run_rg "Monocyte_percentage" \
  "${MUNGED_DIR}/Monocyte_percentage_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.monocyte_percentage.munge.sumstats.gz"

run_rg "Neutrophil_percentage" \
  "${MUNGED_DIR}/Neutrophil_percentage_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.neutrophill_percentage.munge.sumstats.gz"

run_rg "Weight" \
  "${MUNGED_DIR}/Weight_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.weight.munge.sumstats.gz"


# Binary traits (need prevalences)
run_rg "Asthma" \
  "${MUNGED_DIR}/Asthma_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.asthma_doctor_diagnosed.munge.sumstats.gz" \
  --samp-prev 0.10583388,0.127654243 \
  --pop-prev  0.10583388,0.127654243

run_rg "Schizophrenia" \
  "${MUNGED_DIR}/Schizophrenia_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.schizophrenia_ICD10.munge.sumstats.gz" \
  --samp-prev 0.00475155,0.000548182 \
  --pop-prev  0.00475155,0.000548182

run_rg "Type_1_diabetes" \
  "${MUNGED_DIR}/Type_1_diabetes_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.type1_diabetes.munge.sumstats.gz" \
  --samp-prev 0.011693712,0.001614091 \
  --pop-prev  0.011693712,0.001614091

run_rg "Type_2_diabetes" \
  "${MUNGED_DIR}/Type_2_diabetes_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.type2_diabetes.munge.sumstats.gz" \
  --samp-prev 0.115608342,0.002458513 \
  --pop-prev  0.115608342,0.002458513

echo "Completed"

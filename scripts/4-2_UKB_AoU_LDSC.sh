#!/bin/bash
set -euo pipefail

PYTHON="${PYTHON:-python}"
LDSC="/software/ldsc/ldsc.py"

MUNGED_DIR="/munged/gwas_processed"
LDCHR_DIR="/eur_w_ld_chr/"

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


# quantitative traits
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
  "${MUNGED_DIR}/Neale.neutrophil_percentage.munge.sumstats.gz"

run_rg "Weight" \
  "${MUNGED_DIR}/Weight_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.weight.munge.sumstats.gz"

run_rg "White_blood_cell_count" \
  "${MUNGED_DIR}/White_blood_cell_count_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.white_blood_cell_count.munge.sumstats.gz"

run_rg "Red_blood_cell_count" \
  "${MUNGED_DIR}/Red_blood_cell_count_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.red_blood_cell_count.munge.sumstats.gz"

run_rg "Mean_corpuscular_hemoglobin" \
  "${MUNGED_DIR}/Mean_corpuscular_hemoglobin_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.mean_corpuscular_hemoglobin.munge.sumstats.gz"

# binary traits
run_rg "Alzheimers" \
  "${MUNGED_DIR}/Alzheimers_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.alzheimers.munge.sumstats.gz" \
  --samp-prev 0.0022742201,0.0003294628 \
  --pop-prev  0.0022742201,0.0003294628

run_rg "Asthma" \
  "${MUNGED_DIR}/Asthma_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.asthma.munge.sumstats.gz" \
  --samp-prev 0.1286283732,0.0046872318 \
  --pop-prev  0.1286283732,0.0046872318

run_rg "Schizophrenia" \
  "${MUNGED_DIR}/Schizophrenia_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.schizophrenia_ICD10.munge.sumstats.gz" \
  --samp-prev 0.0052631906,0.0005481819 \
  --pop-prev  0.0052631906,0.0005481819

run_rg "Type_1_diabetes" \
  "${MUNGED_DIR}/Type_1_diabetes_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.type1_diabetes.munge.sumstats.gz" \
  --samp-prev 0.0141059285,0.0016140910 \
  --pop-prev  0.0141059285,0.0016140910

run_rg "Type_2_diabetes" \
  "${MUNGED_DIR}/Type_2_diabetes_AoU.munge.sumstats.gz" \
  "${MUNGED_DIR}/Neale.type2_diabetes.munge.sumstats.gz" \
  --samp-prev 0.1454048126,0.0024585126 \
  --pop-prev  0.1454048126,0.0024585126

#!/bin/bash
### Run ASH on all GWAS sumstats

RSCRIPT="${RSCRIPT:-Rscript}"
ASH_SCRIPT="/scripts/7-1_Run_ash_GWAS.R"

# Output 
UKB_BASE="/UKB_sumstats/"
AOU_BASE="/AoU_sumstats/"
FG_BASE="/FG_sumstats"

UKB_IN="${UKB_BASE}/cleaned"
AOU_IN="${AOU_BASE}/cleaned"
FG_IN="${FG_BASE}/cleaned"

UKB_OUT="${UKB_BASE}/ash_out"
AOU_OUT="${AOU_BASE}/ash_out"
FG_OUT="${FG_BASE}/ash_out"

mkdir -p "${UKB_OUT}/logs" "${AOU_OUT}/logs" "${FG_OUT}/logs"

run_ash () {
  local in_file="$1"
  local out_dir="$2"
  local log_file="$3"

  echo "ASH: $(basename "$in_file")"
  "$RSCRIPT" "$ASH_SCRIPT" -f "$in_file" -o "$out_dir" >> "$log_file" 2>&1
}


# UKB 
run_ash "${UKB_IN}/Neale.basophil_percentage.cleaned.txt"     "$UKB_OUT" "${UKB_OUT}/logs/Neale.basophil_percentage.log"
run_ash "${UKB_IN}/Neale.BMI.cleaned.txt"                     "$UKB_OUT" "${UKB_OUT}/logs/Neale.BMI.log"
run_ash "${UKB_IN}/Neale.monocyte_percentage.cleaned.txt"     "$UKB_OUT" "${UKB_OUT}/logs/Neale.monocyte_percentage.log"
run_ash "${UKB_IN}/Neale.neutrophill_percentage.cleaned.txt"  "$UKB_OUT" "${UKB_OUT}/logs/Neale.neutrophill_percentage.log"
run_ash "${UKB_IN}/Neale.schizophrenia_ICD10.cleaned.txt"     "$UKB_OUT" "${UKB_OUT}/logs/Neale.schizophrenia_ICD10.log"
run_ash "${UKB_IN}/Neale.standing_height.cleaned.txt"         "$UKB_OUT" "${UKB_OUT}/logs/Neale.standing_height.log"
run_ash "${UKB_IN}/Neale.type1_diabetes.cleaned.txt"          "$UKB_OUT" "${UKB_OUT}/logs/Neale.type1_diabetes.log"
run_ash "${UKB_IN}/Neale.type2_diabetes.cleaned.txt"          "$UKB_OUT" "${UKB_OUT}/logs/Neale.type2_diabetes.log"
run_ash "${UKB_IN}/Neale.weight.cleaned.txt"                  "$UKB_OUT" "${UKB_OUT}/logs/Neale.weight.log"

#AoU
run_ash "${AOU_IN}/Asthma_AoU.cleaned.txt"                "$AOU_OUT" "${AOU_OUT}/logs/Asthma_AoU.log"
run_ash "${AOU_IN}/Basophil_percentage_AoU.cleaned.txt"   "$AOU_OUT" "${AOU_OUT}/logs/Basophil_percentage_AoU.log"
run_ash "${AOU_IN}/BMI_AoU.cleaned.txt"                   "$AOU_OUT" "${AOU_OUT}/logs/BMI_AoU.log"
run_ash "${AOU_IN}/Height_AoU.cleaned.txt"                "$AOU_OUT" "${AOU_OUT}/logs/Height_AoU.log"
run_ash "${AOU_IN}/Monocyte_percentage_AoU.cleaned.txt"   "$AOU_OUT" "${AOU_OUT}/logs/Monocyte_percentage_AoU.log"
run_ash "${AOU_IN}/Neutrophil_percentage_AoU.cleaned.txt" "$AOU_OUT" "${AOU_OUT}/logs/Neutrophil_percentage_AoU.log"
run_ash "${AOU_IN}/Schizophrenia_AoU.cleaned.txt"         "$AOU_OUT" "${AOU_OUT}/logs/Schizophrenia_AoU.log"
run_ash "${AOU_IN}/Type_1_diabetes_AoU.cleaned.txt"       "$AOU_OUT" "${AOU_OUT}/logs/Type_1_diabetes_AoU.log"
run_ash "${AOU_IN}/Type_2_diabetes_AoU.cleaned.txt"       "$AOU_OUT" "${AOU_OUT}/logs/Type_2_diabetes_AoU.log"
run_ash "${AOU_IN}/Weight_AoU.cleaned.txt"                "$AOU_OUT" "${AOU_OUT}/logs/Weight_AoU.log"

#FG
run_ash "${FG_IN}/finngen_R12_F5_SCHZPHR.cleaned.txt"         "$FG_OUT" "${FG_OUT}/logs/finngen_R12_F5_SCHZPHR.log"
run_ash "${FG_IN}/finngen_R12_G6_ALZHEIMER.cleaned.txt"       "$FG_OUT" "${FG_OUT}/logs/finngen_R12_G6_ALZHEIMER.log"
run_ash "${FG_IN}/finngen_R12_J10_ASTHMA_EXMORE.cleaned.txt"  "$FG_OUT" "${FG_OUT}/logs/finngen_R12_J10_ASTHMA_EXMORE.log"
run_ash "${FG_IN}/finngen_R12_T1D.cleaned.txt"                "$FG_OUT" "${FG_OUT}/logs/finngen_R12_T1D.log"
run_ash "${FG_IN}/finngen_R12_T2D.cleaned.txt"                "$FG_OUT" "${FG_OUT}/logs/finngen_R12_T2D.log"

echo "Done."

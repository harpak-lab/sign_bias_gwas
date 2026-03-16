#!/bin/bash

#set DIRs
PYTHON="${PYTHON:-python}"
RSCRIPT="${RSCRIPT:-Rscript}"

RAW_DIR="/GWAS/UKB_groups/raw"
COMBINED_DIR="/GWAS/UKB_groups/combined"

PREP_SCRIPT="/scripts/4-3_Prep_UKB_group_sumstats.R"
MUNGE="/software/ldsc/munge_sumstats.py"
LDSC="/software/ldsc/ldsc.py"

OUT_DIR="${OUT_DIR:-/munged/gwas_processed}"
RG_DIR="${RG_DIR:-ldsc_rg_out}"

MERGE_ALLELES="/eur_w_ld_chr/w_hm3.snplist"
LDCHR_DIR="/eur_w_ld_chr/"
CHUNKSIZE="${CHUNKSIZE:-500000}"


combine_linear () {
  local group="$1"
  local trait="$2"
  local outfile="$3"

  mapfile -t chr_files < <(
    find "$RAW_DIR" -maxdepth 1 -type f \
      -name "gwas_${group}_chr*.${trait}.glm.linear" | sort -V
  )

  local tmp="${outfile%.gz}"
  awk 'FNR==1 && NR!=1 {next} {print}' "${chr_files[@]}" > "$tmp"
  gzip -f "$tmp"
}

run_prep () {
  local infile="$1"
  local outprefix="$2"

  echo "preprocessing: $(basename "$infile") -> ${outprefix}_filtered.txt"

  "$RSCRIPT" "$PREP_SCRIPT" \
    --file "$infile" \
    --out "$outprefix"
}

run_munge () {
  local sumstats="$1"
  local outprefix="$2"
  local a1="$3"
  local a2="$4"
  shift 4

  echo "munging: $(basename "$sumstats") -> ${outprefix}.sumstats.gz"

  "$PYTHON" "$MUNGE" \
    --sumstats "$sumstats" \
    --a1 "$a1" \
    --a2 "$a2" \
    --snp rsid \
    --p pval \
    --merge-alleles "$MERGE_ALLELES" \
    --chunksize "$CHUNKSIZE" \
    --out "$outprefix" \
    "$@"
}

run_rg () {
  local trait="$1"
  local f1="$2"
  local f2="$3"
  shift 3

  echo "LDSC: $trait"

  "$PYTHON" "$LDSC" \
    --rg "${f1},${f2}" \
    --ref-ld-chr "$LDCHR_DIR" \
    --w-ld-chr "$LDCHR_DIR" \
    --out "${RG_DIR}/${trait}" \
    "$@"
}

#allele cols
UKB_A1="a1"
UKB_A2="a2"


#combine

combine_linear group1 height      "$COMBINED_DIR/UKB_group1_height.tsv.gz"
combine_linear group2 height      "$COMBINED_DIR/UKB_group2_height.tsv.gz"

combine_linear group1 weight      "$COMBINED_DIR/UKB_group1_weight.tsv.gz"
combine_linear group2 weight      "$COMBINED_DIR/UKB_group2_weight.tsv.gz"

combine_linear group1 bmi         "$COMBINED_DIR/UKB_group1_BMI.tsv.gz"
combine_linear group2 bmi         "$COMBINED_DIR/UKB_group2_BMI.tsv.gz"

combine_linear group1 mono_pct    "$COMBINED_DIR/UKB_group1_mono_pct.tsv.gz"
combine_linear group2 mono_pct    "$COMBINED_DIR/UKB_group2_mono_pct.tsv.gz"

combine_linear group1 baso_pct    "$COMBINED_DIR/UKB_group1_baso_pct.tsv.gz"
combine_linear group2 baso_pct    "$COMBINED_DIR/UKB_group2_baso_pct.tsv.gz"

combine_linear group1 neutro_pct  "$COMBINED_DIR/UKB_group1_neutro_pct.tsv.gz"
combine_linear group2 neutro_pct  "$COMBINED_DIR/UKB_group2_neutro_pct.tsv.gz"

combine_linear group1 wbc         "$COMBINED_DIR/UKB_group1_wbc.tsv.gz"
combine_linear group2 wbc         "$COMBINED_DIR/UKB_group2_wbc.tsv.gz"

combine_linear group1 rbc         "$COMBINED_DIR/UKB_group1_rbc.tsv.gz"
combine_linear group2 rbc         "$COMBINED_DIR/UKB_group2_rbc.tsv.gz"

combine_linear group1 mch         "$COMBINED_DIR/UKB_group1_mch.tsv.gz"
combine_linear group2 mch         "$COMBINED_DIR/UKB_group2_mch.tsv.gz"

combine_linear group1 AD          "$COMBINED_DIR/UKB_group1_AD.tsv.gz"
combine_linear group2 AD          "$COMBINED_DIR/UKB_group2_AD.tsv.gz"

combine_linear group1 ASTHMA      "$COMBINED_DIR/UKB_group1_ASTHMA.tsv.gz"
combine_linear group2 ASTHMA      "$COMBINED_DIR/UKB_group2_ASTHMA.tsv.gz"

combine_linear group1 SCZ         "$COMBINED_DIR/UKB_group1_SCZ.tsv.gz"
combine_linear group2 SCZ         "$COMBINED_DIR/UKB_group2_SCZ.tsv.gz"

combine_linear group1 T1D         "$COMBINED_DIR/UKB_group1_T1D.tsv.gz"
combine_linear group2 T1D         "$COMBINED_DIR/UKB_group2_T1D.tsv.gz"

combine_linear group1 T2D         "$COMBINED_DIR/UKB_group1_T2D.tsv.gz"
combine_linear group2 T2D         "$COMBINED_DIR/UKB_group2_T2D.tsv.gz"


# prep (pre-munge)

run_prep "$COMBINED_DIR/UKB_group1_height.tsv.gz"      "$OUT_DIR/UKB_group1_height"
run_prep "$COMBINED_DIR/UKB_group2_height.tsv.gz"      "$OUT_DIR/UKB_group2_height"

run_prep "$COMBINED_DIR/UKB_group1_weight.tsv.gz"      "$OUT_DIR/UKB_group1_weight"
run_prep "$COMBINED_DIR/UKB_group2_weight.tsv.gz"      "$OUT_DIR/UKB_group2_weight"

run_prep "$COMBINED_DIR/UKB_group1_BMI.tsv.gz"         "$OUT_DIR/UKB_group1_BMI"
run_prep "$COMBINED_DIR/UKB_group2_BMI.tsv.gz"         "$OUT_DIR/UKB_group2_BMI"

run_prep "$COMBINED_DIR/UKB_group1_mono_pct.tsv.gz"   "$OUT_DIR/UKB_group1_mono_pct"
run_prep "$COMBINED_DIR/UKB_group2_mono_pct.tsv.gz"   "$OUT_DIR/UKB_group2_mono_pct"

run_prep "$COMBINED_DIR/UKB_group1_baso_pct.tsv.gz"   "$OUT_DIR/UKB_group1_baso_pct"
run_prep "$COMBINED_DIR/UKB_group2_baso_pct.tsv.gz"   "$OUT_DIR/UKB_group2_baso_pct"

run_prep "$COMBINED_DIR/UKB_group1_neutro_pct.tsv.gz" "$OUT_DIR/UKB_group1_neutro_pct"
run_prep "$COMBINED_DIR/UKB_group2_neutro_pct.tsv.gz" "$OUT_DIR/UKB_group2_neutro_pct"

run_prep "$COMBINED_DIR/UKB_group1_wbc.tsv.gz"         "$OUT_DIR/UKB_group1_wbc"
run_prep "$COMBINED_DIR/UKB_group2_wbc.tsv.gz"         "$OUT_DIR/UKB_group2_wbc"

run_prep "$COMBINED_DIR/UKB_group1_rbc.tsv.gz"         "$OUT_DIR/UKB_group1_rbc"
run_prep "$COMBINED_DIR/UKB_group2_rbc.tsv.gz"         "$OUT_DIR/UKB_group2_rbc"

run_prep "$COMBINED_DIR/UKB_group1_mch.tsv.gz"         "$OUT_DIR/UKB_group1_mch"
run_prep "$COMBINED_DIR/UKB_group2_mch.tsv.gz"         "$OUT_DIR/UKB_group2_mch"

run_prep "$COMBINED_DIR/UKB_group1_AD.tsv.gz"          "$OUT_DIR/UKB_group1_AD"
run_prep "$COMBINED_DIR/UKB_group2_AD.tsv.gz"          "$OUT_DIR/UKB_group2_AD"

run_prep "$COMBINED_DIR/UKB_group1_ASTHMA.tsv.gz"      "$OUT_DIR/UKB_group1_ASTHMA"
run_prep "$COMBINED_DIR/UKB_group2_ASTHMA.tsv.gz"      "$OUT_DIR/UKB_group2_ASTHMA"

run_prep "$COMBINED_DIR/UKB_group1_SCZ.tsv.gz"         "$OUT_DIR/UKB_group1_SCZ"
run_prep "$COMBINED_DIR/UKB_group2_SCZ.tsv.gz"         "$OUT_DIR/UKB_group2_SCZ"

run_prep "$COMBINED_DIR/UKB_group1_T1D.tsv.gz"         "$OUT_DIR/UKB_group1_T1D"
run_prep "$COMBINED_DIR/UKB_group2_T1D.tsv.gz"         "$OUT_DIR/UKB_group2_T1D"

run_prep "$COMBINED_DIR/UKB_group1_T2D.tsv.gz"         "$OUT_DIR/UKB_group1_T2D"
run_prep "$COMBINED_DIR/UKB_group2_T2D.tsv.gz"         "$OUT_DIR/UKB_group2_T2D"

#munge step (LDSC)

run_munge "$OUT_DIR/UKB_group1_height_filtered.txt" "$OUT_DIR/UKB_group1_height.munge" "$UKB_A1" "$UKB_A2" --N 180079
run_munge "$OUT_DIR/UKB_group2_height_filtered.txt" "$OUT_DIR/UKB_group2_height.munge" "$UKB_A1" "$UKB_A2" --N 180067

run_munge "$OUT_DIR/UKB_group1_weight_filtered.txt" "$OUT_DIR/UKB_group1_weight.munge" "$UKB_A1" "$UKB_A2" --N 179944
run_munge "$OUT_DIR/UKB_group2_weight_filtered.txt" "$OUT_DIR/UKB_group2_weight.munge" "$UKB_A1" "$UKB_A2" --N 179929

run_munge "$OUT_DIR/UKB_group1_BMI_filtered.txt" "$OUT_DIR/UKB_group1_BMI.munge" "$UKB_A1" "$UKB_A2" --N 179880
run_munge "$OUT_DIR/UKB_group2_BMI_filtered.txt" "$OUT_DIR/UKB_group2_BMI.munge" "$UKB_A1" "$UKB_A2" --N 179861

run_munge "$OUT_DIR/UKB_group1_mono_pct_filtered.txt" "$OUT_DIR/UKB_group1_mono_pct.munge" "$UKB_A1" "$UKB_A2" --N 174890
run_munge "$OUT_DIR/UKB_group2_mono_pct_filtered.txt" "$OUT_DIR/UKB_group2_mono_pct.munge" "$UKB_A1" "$UKB_A2" --N 174739

run_munge "$OUT_DIR/UKB_group1_baso_pct_filtered.txt" "$OUT_DIR/UKB_group1_baso_pct.munge" "$UKB_A1" "$UKB_A2" --N 174890
run_munge "$OUT_DIR/UKB_group2_baso_pct_filtered.txt" "$OUT_DIR/UKB_group2_baso_pct.munge" "$UKB_A1" "$UKB_A2" --N 174739

run_munge "$OUT_DIR/UKB_group1_neutro_pct_filtered.txt" "$OUT_DIR/UKB_group1_neutro_pct.munge" "$UKB_A1" "$UKB_A2" --N 174890
run_munge "$OUT_DIR/UKB_group2_neutro_pct_filtered.txt" "$OUT_DIR/UKB_group2_neutro_pct.munge" "$UKB_A1" "$UKB_A2" --N 174739

run_munge "$OUT_DIR/UKB_group1_wbc_filtered.txt" "$OUT_DIR/UKB_group1_wbc.munge" "$UKB_A1" "$UKB_A2" --N 175208
run_munge "$OUT_DIR/UKB_group2_wbc_filtered.txt" "$OUT_DIR/UKB_group2_wbc.munge" "$UKB_A1" "$UKB_A2" --N 175029

run_munge "$OUT_DIR/UKB_group1_rbc_filtered.txt" "$OUT_DIR/UKB_group1_rbc.munge" "$UKB_A1" "$UKB_A2" --N 175211
run_munge "$OUT_DIR/UKB_group2_rbc_filtered.txt" "$OUT_DIR/UKB_group2_rbc.munge" "$UKB_A1" "$UKB_A2" --N 175031

run_munge "$OUT_DIR/UKB_group1_mch_filtered.txt" "$OUT_DIR/UKB_group1_mch.munge" "$UKB_A1" "$UKB_A2" --N 175208
run_munge "$OUT_DIR/UKB_group2_mch_filtered.txt" "$OUT_DIR/UKB_group2_mch.munge" "$UKB_A1" "$UKB_A2" --N 175031

# Binary
run_munge "$OUT_DIR/UKB_group1_AD_filtered.txt" "$OUT_DIR/UKB_group1_AD.munge" "$UKB_A1" "$UKB_A2" --N-cas 971 --N-con 180232
run_munge "$OUT_DIR/UKB_group2_AD_filtered.txt" "$OUT_DIR/UKB_group2_AD.munge" "$UKB_A1" "$UKB_A2" --N-cas 903 --N-con 180250

run_munge "$OUT_DIR/UKB_group1_ASTHMA_filtered.txt" "$OUT_DIR/UKB_group1_ASTHMA.munge" "$UKB_A1" "$UKB_A2" --N-cas 5438 --N-con 179105
run_munge "$OUT_DIR/UKB_group2_ASTHMA_filtered.txt" "$OUT_DIR/UKB_group2_ASTHMA.munge" "$UKB_A1" "$UKB_A2" --N-cas 5501 --N-con 179090

run_munge "$OUT_DIR/UKB_group1_SCZ_filtered.txt" "$OUT_DIR/UKB_group1_SCZ.munge" "$UKB_A1" "$UKB_A2" --N-cas 572 --N-con 180332
run_munge "$OUT_DIR/UKB_group2_SCZ_filtered.txt" "$OUT_DIR/UKB_group2_SCZ.munge" "$UKB_A1" "$UKB_A2" --N-cas 619 --N-con 180321

run_munge "$OUT_DIR/UKB_group1_T1D_filtered.txt" "$OUT_DIR/UKB_group1_T1D.munge" "$UKB_A1" "$UKB_A2" --N-cas 1477 --N-con 180105
run_munge "$OUT_DIR/UKB_group2_T1D_filtered.txt" "$OUT_DIR/UKB_group2_T1D.munge" "$UKB_A1" "$UKB_A2" --N-cas 1433 --N-con 180117

run_munge "$OUT_DIR/UKB_group1_T2D_filtered.txt" "$OUT_DIR/UKB_group1_T2D.munge" "$UKB_A1" "$UKB_A2" --N-cas 2777 --N-con 179778
run_munge "$OUT_DIR/UKB_group2_T2D_filtered.txt" "$OUT_DIR/UKB_group2_T2D.munge" "$UKB_A1" "$UKB_A2" --N-cas 2749 --N-con 179786


#LDSC rg 

#ukb 1 vs ukb 2
run_rg "Height_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_height.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_height.munge.sumstats.gz"

run_rg "Weight_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_weight.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_weight.munge.sumstats.gz"

run_rg "BMI_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_BMI.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_BMI.munge.sumstats.gz"

run_rg "Monocyte_percentage_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_monocyte_percentage.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_monocyte_percentage.munge.sumstats.gz"

run_rg "Basophil_percentage_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_basophil_percentage.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_basophil_percentage.munge.sumstats.gz"

run_rg "Neutrophil_percentage_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_neutrophil_percentage.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_neutrophil_percentage.munge.sumstats.gz"

run_rg "White_blood_cell_count_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_wbc.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_wbc.munge.sumstats.gz"

run_rg "Red_blood_cell_count_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_rbc.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_rbc.munge.sumstats.gz"

run_rg "Mean_corpuscular_hemoglobin_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_mch.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_mch.munge.sumstats.gz"

# Quantitative: UKB1 vs AoU
run_rg "Height_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_height.munge.sumstats.gz" \
  "$OUT_DIR/Height_AoU.munge.sumstats.gz"

run_rg "Weight_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_weight.munge.sumstats.gz" \
  "$OUT_DIR/Weight_AoU.munge.sumstats.gz"

run_rg "BMI_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_BMI.munge.sumstats.gz" \
  "$OUT_DIR/BMI_AoU.munge.sumstats.gz"

run_rg "Monocyte_percentage_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_monocyte_percentage.munge.sumstats.gz" \
  "$OUT_DIR/Monocyte_percentage_AoU.munge.sumstats.gz"

run_rg "Basophil_percentage_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_basophil_percentage.munge.sumstats.gz" \
  "$OUT_DIR/Basophil_percentage_AoU.munge.sumstats.gz"

run_rg "Neutrophil_percentage_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_neutrophil_percentage.munge.sumstats.gz" \
  "$OUT_DIR/Neutrophil_percentage_AoU.munge.sumstats.gz"

run_rg "White_blood_cell_count_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_wbc.munge.sumstats.gz" \
  "$OUT_DIR/White_blood_cell_count_AoU.munge.sumstats.gz"

run_rg "Red_blood_cell_count_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_rbc.munge.sumstats.gz" \
  "$OUT_DIR/Red_blood_cell_count_AoU.munge.sumstats.gz"

run_rg "Mean_corpuscular_hemoglobin_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_mch.munge.sumstats.gz" \
  "$OUT_DIR/Mean_corpuscular_hemoglobin_AoU.munge.sumstats.gz"

# Quantitative: UKB2 vs AoU
run_rg "Height_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_height.munge.sumstats.gz" \
  "$OUT_DIR/Height_AoU.munge.sumstats.gz"

run_rg "Weight_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_weight.munge.sumstats.gz" \
  "$OUT_DIR/Weight_AoU.munge.sumstats.gz"

run_rg "BMI_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_BMI.munge.sumstats.gz" \
  "$OUT_DIR/BMI_AoU.munge.sumstats.gz"

run_rg "Monocyte_percentage_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_monocyte_percentage.munge.sumstats.gz" \
  "$OUT_DIR/Monocyte_percentage_AoU.munge.sumstats.gz"

run_rg "Basophil_percentage_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_basophil_percentage.munge.sumstats.gz" \
  "$OUT_DIR/Basophil_percentage_AoU.munge.sumstats.gz"

run_rg "Neutrophil_percentage_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_neutrophil_percentage.munge.sumstats.gz" \
  "$OUT_DIR/Neutrophil_percentage_AoU.munge.sumstats.gz"

run_rg "White_blood_cell_count_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_wbc.munge.sumstats.gz" \
  "$OUT_DIR/White_blood_cell_count_AoU.munge.sumstats.gz"

run_rg "Red_blood_cell_count_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_rbc.munge.sumstats.gz" \
  "$OUT_DIR/Red_blood_cell_count_AoU.munge.sumstats.gz"

run_rg "Mean_corpuscular_hemoglobin_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_mch.munge.sumstats.gz" \
  "$OUT_DIR/Mean_corpuscular_hemoglobin_AoU.munge.sumstats.gz"

# Binary: UKB1 vs UKB2
run_rg "Alzheimers_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_AD.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_AD.munge.sumstats.gz" \
  --samp-prev 0.0053550518,0.0049853000 \
  --pop-prev 0.0053550518,0.0049853000

run_rg "Asthma_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_ASTHMA.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_ASTHMA.munge.sumstats.gz" \
  --samp-prev 0.0294783867,0.0298049644 \
  --pop-prev 0.0294783867,0.0298049644

run_rg "Schizophrenia_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_SCZ.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_SCZ.munge.sumstats.gz" \
  --samp-prev 0.0031623440,0.0034209456 \
  --pop-prev 0.0031623440,0.0034209456

run_rg "Type_1_diabetes_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_T1D.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_T1D.munge.sumstats.gz" \
  --samp-prev 0.0081353074,0.0078917836 \
  --pop-prev 0.0081353074,0.0078917836

run_rg "Type_2_diabetes_UKB1_vs_UKB2" \
  "$OUT_DIR/UKB_group1_T2D.munge.sumstats.gz" \
  "$OUT_DIR/UKB_group2_T2D.munge.sumstats.gz" \
  --samp-prev 0.0152094681,0.0150571814 \
  --pop-prev 0.0152094681,0.0150571814

# Binary: UKB1 vs AoU
run_rg "Alzheimers_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_AD.munge.sumstats.gz" \
  "$OUT_DIR/Alzheimers_AoU.munge.sumstats.gz" \
  --samp-prev 0.0053550518,0.0022742201 \
  --pop-prev 0.0053550518,0.0022742201

run_rg "Asthma_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_ASTHMA.munge.sumstats.gz" \
  "$OUT_DIR/Asthma_AoU.munge.sumstats.gz" \
  --samp-prev 0.0294783867,0.1286283732 \
  --pop-prev 0.0294783867,0.1286283732

run_rg "Schizophrenia_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_SCZ.munge.sumstats.gz" \
  "$OUT_DIR/Schizophrenia_AoU.munge.sumstats.gz" \
  --samp-prev 0.0031623440,0.0052631906 \
  --pop-prev 0.0031623440,0.0052631906

run_rg "Type_1_diabetes_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_T1D.munge.sumstats.gz" \
  "$OUT_DIR/Type_1_diabetes_AoU.munge.sumstats.gz" \
  --samp-prev 0.0081353074,0.0141059285 \
  --pop-prev 0.0081353074,0.0141059285

run_rg "Type_2_diabetes_UKB1_vs_AoU" \
  "$OUT_DIR/UKB_group1_T2D.munge.sumstats.gz" \
  "$OUT_DIR/Type_2_diabetes_AoU.munge.sumstats.gz" \
  --samp-prev 0.0152094681,0.1454048126 \
  --pop-prev 0.0152094681,0.1454048126

# Binary: UKB2 vs AoU
run_rg "Alzheimers_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_AD.munge.sumstats.gz" \
  "$OUT_DIR/Alzheimers_AoU.munge.sumstats.gz" \
  --samp-prev 0.0049853000,0.0022742201 \
  --pop-prev 0.0049853000,0.0022742201 

run_rg "Asthma_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_ASTHMA.munge.sumstats.gz" \
  "$OUT_DIR/Asthma_AoU.munge.sumstats.gz" \
  --samp-prev 0.0298049644,0.1286283732 \
  --pop-prev 0.0298049644,0.1286283732

run_rg "Schizophrenia_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_SCZ.munge.sumstats.gz" \
  "$OUT_DIR/Schizophrenia_AoU.munge.sumstats.gz" \
  --samp-prev 0.0034209456,0.0052631906 \
  --pop-prev 0.0034209456,0.0052631906

run_rg "Type_1_diabetes_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_T1D.munge.sumstats.gz" \
  "$OUT_DIR/Type_1_diabetes_AoU.munge.sumstats.gz" \
  --samp-prev 0.0078917836,0.0141059285 \
  --pop-prev 0.0078917836,0.0141059285

run_rg "Type_2_diabetes_UKB2_vs_AoU" \
  "$OUT_DIR/UKB_group2_T2D.munge.sumstats.gz" \
  "$OUT_DIR/Type_2_diabetes_AoU.munge.sumstats.gz" \
  --samp-prev 0.0150571814,0.1454048126 \
  --pop-prev 0.0150571814,0.1454048126
#!/bin/bash

#DIRs
PHENO_DIR="${PHENO_DIR:-phenos}"
COVAR_IN="${COVAR_IN:-${PHENO_DIR}/aou_covariates.centered.tsv}"
COVAR_PLINK="${COVAR_PLINK:-${PHENO_DIR}/aou_covariates.centered.plink.tsv}"

GENO_BASE="${GENO_BASE:-plink}"
OUT_BASE="${OUT_BASE:-gwas_out}"
THREADS="${THREADS:-8}"
MAX_PARALLEL_TRAITS="${MAX_PARALLEL_TRAITS:-5}" #run max 5 traits in parallel 

RECODE_BINARY_TO_34="${RECODE_BINARY_TO_34:-1}"

if [[ ! -f "$COVAR_PLINK" ]]; then
  awk -F'\t' -v OFS='\t' '
    NR==1{
      printf "FID\tIID";
      for(i=2;i<=NF;i++){ printf "\t%s", $i }
      printf "\n";
      next
    }
    {
      id=$1;
      printf "%s\t%s", id, id;
      for(i=2;i<=NF;i++){ printf "\t%s", $i }
      printf "\n";
    }
  ' "$COVAR_IN" > "$COVAR_PLINK"
fi

#phenotypes
mapfile -t PHENO_FILES < <(find "$PHENO_DIR" -maxdepth 1 -type f -name "aou_*.plink.tsv" ! -name "*covariates*" | sort)

pick_prefix () {
  local chr="$1"
  local p1="${GENO_BASE}/chr${chr}.filtered"
  local p2="${GENO_BASE}/chr${chr}/chr${chr}.filtered"

  if [[ -f "${p1}.pgen" || -f "${p1}.bed" ]]; then
    echo "$p1"
    return 0
  elif [[ -f "${p2}.pgen" || -f "${p2}.bed" ]]; then
    echo "$p2"
    return 0
  fi
  return 1
}


#recode binrary to run linreg
maybe_recode_pheno_to_34 () {
  local pheno_file="$1"
  local trait="$2"

  if [[ "$RECODE_BINARY_TO_34" -ne 1 ]]; then
    echo "$pheno_file"
    return 0
  fi

  if awk -F'\t' '
      NR==1{next}
      $3=="" || $3=="NA" || $3=="." {next}
      { seen[$3]=1; n++ }
      END{
        if(n==0) exit 1
        for(k in seen){
          if(!(k==0 || k==1 || k==2)) exit 1
        }
        if(seen[0] && seen[2]) exit 1
        exit 0
      }
    ' "$pheno_file"
  then
    local tmpdir="${OUT_BASE}/.tmp_pheno"
    mkdir -p "$tmpdir"
    local tmpfile="${tmpdir}/${trait}.to34.plink.tsv"

    if awk -F'\t' 'NR>1 && $3==0 {found=1; exit} END{exit(found?0:1)}' "$pheno_file"; then
      awk -F'\t' -v OFS='\t' '
        NR==1{print; next}
        {
          if($3==0) $3=3;
          else if($3==1) $3=4;
          print
        }
      ' "$pheno_file" > "$tmpfile"
    else
      awk -F'\t' -v OFS='\t' '
        NR==1{print; next}
        {
          if($3==1) $3=3;
          else if($3==2) $3=4;
          print
        }
      ' "$pheno_file" > "$tmpfile"
    fi

    echo "$tmpfile"
    return 0
  else
    echo "$pheno_file"
    return 0
  fi
}

# map phenotype basename to outnames
merged_trait_name () {
  local trait="$1"
  case "$trait" in
    alzheimers) echo "Alzheimers_AoU" ;;
    asthma) echo "Asthma_AoU" ;;
    basophil_percentage) echo "Basophil_percentage_AoU" ;;
    bmi) echo "BMI_AoU" ;;
    height) echo "Height_AoU" ;;
    monocyte_percentage) echo "Monocyte_percentage_AoU" ;;
    neutrophil_percentage) echo "Neutrophil_percentage_AoU" ;;
    white_blood_cell_count) echo "White_blood_cell_count_AoU" ;;
    red_blood_cell_count) echo "Red_blood_cell_count_AoU" ;;
    mean_corpuscular_hemoglobin) echo "Mean_corpuscular_hemoglobin_AoU" ;;
    schizophrenia) echo "Schizophrenia_AoU" ;;
    t1d) echo "Type_1_diabetes_AoU" ;;
    t2d) echo "Type_2_diabetes_AoU" ;;
    weight) echo "Weight_AoU" ;;
    *)
      echo "${trait}_AoU"
      ;;
  esac
}

mkdir -p "${OUT_BASE}/logs"

run_trait_gwas () {
  local pf="$1"

  local bn trait outdir pheno_use outprefix logfile prefix INFLAG chr
  bn="$(basename "$pf")"
  trait="${bn#aou_}"
  trait="${trait%.plink.tsv}"

  outdir="${OUT_BASE}/${trait}"
  mkdir -p "$outdir"

  pheno_use="$(maybe_recode_pheno_to_34 "$pf" "$trait")"

  for chr in $(seq 1 22); do
    prefix="$(pick_prefix "$chr")" || continue

    if [[ -f "${prefix}.pgen" ]]; then
      INFLAG="--pfile"
    elif [[ -f "${prefix}.bed" ]]; then
      INFLAG="--bfile"
    else
      continue
    fi

    outprefix="${outdir}/chr${chr}"
    logfile="${OUT_BASE}/logs/${trait}.chr${chr}.log"

    plink2 \
      ${INFLAG} "${prefix}" \
      --pheno "${pheno_use}" \
      --covar "${COVAR_PLINK}" \
      --covar-variance-standardize \
      --glm hide-covar omit-ref \
      --vif 80 \
      --threads "${THREADS}" \
      --out "${outprefix}" \
      > "${logfile}" 2>&1
  done

  trait_dir="${OUT_BASE}/${trait}"

  mapfile -t chr_files < <(
    find "${trait_dir}" -maxdepth 1 -type f -name "chr*.glm.linear" | sort -V
  )

  merged_name="$(merged_trait_name "$trait")"
  tmpfile="${OUT_BASE}/${merged_name}.tsv"

  awk 'FNR==1 && NR!=1 {next} {print}' "${chr_files[@]}" > "${tmpfile}"
  gzip -f "${tmpfile}"
}

for pf in "${PHENO_FILES[@]}"; do
  while [ "$(jobs -r | wc -l)" -ge "$MAX_PARALLEL_TRAITS" ]; do
    sleep 2
  done

  run_trait_gwas "$pf" &
done

wait
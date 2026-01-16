#!/bin/bash
## Collect phenotype files, covariates, and run GWAS (plink 2.0)
########

## Set up directories
PHENO_DIR="${PHENO_DIR:-phenos}"
COVAR_IN="${COVAR_IN:-${PHENO_DIR}/aou_covariates.centered.tsv}"
COVAR_PLINK="${COVAR_PLINK:-${PHENO_DIR}/aou_covariates.centered.plink.tsv}"

GENO_BASE="${GENO_BASE:-plink}"
OUT_BASE="${OUT_BASE:-gwas_out}"
THREADS="${THREADS:-8}"

# Recode binary phenotypes (so --glm runs linear regression)
RECODE_BINARY_TO_34="${RECODE_BINARY_TO_34:-1}"

# Format covariate file 
if [[ ! -f "$COVAR_PLINK" ]]; then
  echo "Creating PLINK covariate file: $COVAR_PLINK"
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

# Phenotype files
mapfile -t PHENO_FILES < <(find "$PHENO_DIR" -maxdepth 1 -type f -name "aou_*.plink.tsv" ! -name "*covariates*" | sort)
if [[ "${#PHENO_FILES[@]}" -eq 0 ]]; then
  echo "no matching ${PHENO_DIR}/aou_*.plink.tsv" >&2
  exit 2
fi

# Chromosome prefix
pick_prefix () {
  local chr="$1"
  local p1="${GENO_BASE}/chr${chr}.filtered"
  local p2="${GENO_BASE}/chr${chr}/chr${chr}.filtered"
  return 1
}

# Recode binary to run linreg
maybe_recode_pheno_to_34 () {
  local pheno_file="$1"
  local trait="$2"

  if [[ "$RECODE_BINARY_TO_34" -ne 1 ]]; then
    echo "$pheno_file"; return 0
  fi

  # Detect if all non-missing PHENO values
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

    echo "$tmpfile"; return 0
  else
    echo "$pheno_file"; return 0
  fi
}

for chr in $(seq 1 22); do
  if ! prefix="$(pick_prefix "$chr")"; then
    echo "WARN: (skipping chr${chr})"
    continue
  fi

  if [[ -f "${prefix}.pgen" ]]; then
    INFLAG="--pfile"
  elif [[ -f "${prefix}.bed" ]]; then
    INFLAG="--bfile"
  else
    echo "(skipping)"
    continue
  fi

  for pf in "${PHENO_FILES[@]}"; do
    bn="$(basename "$pf")"
    trait="${bn#aou_}"
    trait="${trait%.plink.tsv}"

    outdir="${OUT_BASE}/${trait}"
    mkdir -p "$outdir"

    pheno_use="$(maybe_recode_pheno_to_34 "$pf" "$trait")"
    outprefix="${outdir}/chr${chr}"
    logfile="${OUT_BASE}/logs/${trait}.chr${chr}.log"

    echo "chr${chr} trait=${trait} input=${prefix} pheno=$(basename "$pheno_use")"

    plink2 \
      ${INFLAG} "${prefix}" \
      --pheno "${pheno_use}" \
      --covar "${COVAR_PLINK}" \
      --covar-variance-standardize \
      --glm hide-covar omit-ref \
      --threads "${THREADS}" \
      --out "${outprefix}" \
      > "${logfile}" 2>&1
  done
done

echo
echo "Done."

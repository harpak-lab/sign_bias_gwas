#!/bin/bash

# QC AoU ACAF variants
######################

PGEN_BASE_DIR="${1:-./plink}"
KEEP_FILE="${2:-./wb_samples_filtered.txt}" #ancestry matched samples

# Chromosomes to run
CHROMS=()
for c in $(seq 1 22); do CHROMS+=("$c"); done

echo "PGEN_BASE_DIR = $PGEN_BASE_DIR"
echo "KEEP_FILE     = $KEEP_FILE"
echo

for chr in "${CHROMS[@]}"; do
  chr_dir="${PGEN_BASE_DIR}/chr${chr}"
  prefix="${chr_dir}/chr${chr}"

  # Find a .pgen 
  if [[ ! -f "${prefix}.pgen" ]]; then
    echo "missing ${prefix}.pgen (skipping chr${chr})"
    continue
  fi

  echo "Running on chr${chr}..."
  (
    cd "$chr_dir"
    plink2 \
      --pfile "chr${chr}" \
      --keep "$(realpath "$KEEP_FILE")" \
      --hwe 1e-12 midp keep-fewhet \
      --geno 0.05 \
      --make-bed \
      --out "chr${chr}.filtered" \
      --allow-extra-chr
  )
done

echo
echo "Done"

#!/bin/bash
### Run all binning and MAF bin alculations

RSCRIPT="${RSCRIPT:-Rscript}"
SCRIPT="7-3_Bin_and_Calc_sign_bias.R"

# Where ASH outputs live (input files for bin/threshold)
UKB_ASH="/UKB/ash_out"
AOU_ASH="/AoU/ash_out"
FG_ASH="/FG/ash_out"

# Output
OUT_BASE="/binned_out/"
mkdir -p "$OUT_BASE"/{UKB,AOU,FG}/{bin,threshold,logs}

TYPE="MAF"
METHODS=("sig" "random")
THRESHOLDS=(0.001 0.01 0.1)

run_one () {
  local biobank="$1"
  local mode="$2"         # bin or threshold
  local method="$3"       # sig or random
  local infile="$4"
  local outdir="$5"
  local logfile="$6"
  shift 6

  mkdir -p "$outdir"
  (
    cd "$outdir"
    echo "==> $(date) biobank=${biobank} mode=${mode} method=${method} file=$(basename "$infile")"
    "$RSCRIPT" "$SCRIPT" \
      --file "$infile" \
      --type "$TYPE" \
      --method "$method" \
      --mode "$mode" \
      "$@"
  ) >> "$logfile" 2>&1
}

collect_files () {
  local dir="$1"
  find "$dir" -type f \( -name "*.txt" -o -name "*.tsv" -o -name "*.txt.gz" -o -name "*.tsv.gz" \) \
    ! -path "*/logs/*" \
    | sort
}

run_biobank () {
  local biobank="$1"
  local indir="$2"
  local out_bin="${OUT_BASE}/${biobank}/bin"
  local out_thr="${OUT_BASE}/${biobank}/threshold"
  local logdir="${OUT_BASE}/${biobank}/logs"

  mapfile -t files < <(collect_files "$indir")
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "WARN: no ASH output files found under $indir for $biobank" >&2
    return 0
  fi

  for f in "${files[@]}"; do
    base="$(basename "$f")"
    base_noext="${base%.*}"
    base_noext="${base_noext%.txt}"
    base_noext="${base_noext%.tsv}"

    for method in "${METHODS[@]}"; do
      # 1) BIN MODE
      run_one "$biobank" "bin" "$method" "$f" "$out_bin" "${logdir}/${base_noext}.${method}.bin.log"

      # 2) THRESHOLD MODE at MAF < T
      for T in "${THRESHOLDS[@]}"; do
        run_one "$biobank" "threshold" "$method" "$f" "$out_thr" "${logdir}/${base_noext}.${method}.thr_${T}.log" --low "$T"
      done
    done
  done
}

echo "Starting"
echo "TYPE=$TYPE"
echo "METHODS: ${METHODS[*]}"
echo "Thresholds: ${THRESHOLDS[*]}"
echo

run_biobank "UKB" "$UKB_ASH"
run_biobank "AOU" "$AOU_ASH"
run_biobank "FG"  "$FG_ASH"

echo
echo "Done."

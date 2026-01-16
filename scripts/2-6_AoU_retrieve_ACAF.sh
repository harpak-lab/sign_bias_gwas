#!/bin/bash
## Get AOU data (ACAF)
######################

#Variant data

CDR_STORAGE_PATH="${CDR_STORAGE_PATH:-gs://fc-aou-datasets-controlled/v8}"

SRC="${CDR_STORAGE_PATH}/wgs/short_read/snpindel/acaf_threshold/pgen"
DEST_BASE="${1:-${PWD}/plink}"
mkdir -p "$DEST_BASE"

# Chromosomes to fetch
CHROMS=()
for c in $(seq 1 22); do CHROMS+=("$c"); done

gcs_dir_exists () {
  local d="$1"
  gsutil -u "$GOOGLE_PROJECT" ls "${d}/" >/dev/null 2>&1
}

if gcs_dir_exists "${SRC}/chr1"; then
  echo "Detected per-chromosome subdirectories under $SRC (e.g., $SRC/chr1/)."
  for chr in "${CHROMS[@]}"; do
    echo "Copying chr${chr} directory -> ${DEST_BASE}/chr${chr}"
    mkdir -p "${DEST_BASE}/chr${chr}"
    gsutil -m -u "$GOOGLE_PROJECT" cp -r "${SRC}/chr${chr}/"* "${DEST_BASE}/chr${chr}/"
  done
else
  tmp_all="$(mktemp)"
  tmp_pfiles="$(mktemp)"
  trap 'rm -f "$tmp_all" "$tmp_pfiles"' EXIT

  # List once 
  gsutil -u "$GOOGLE_PROJECT" ls -r "${SRC}/**" > "$tmp_all"

  # Keep pgen/pvar/psam 
  grep -E '\.(pgen|pvar|psam)(\..+)?$' "$tmp_all" > "$tmp_pfiles" || true

  if [[ ! -s "$tmp_pfiles" ]]; then
    echo "ERROR: No .pgen/.pvar/.psam files found under: $SRC"
    exit 2
  fi

  for chr in "${CHROMS[@]}"; do
    mkdir -p "${DEST_BASE}/chr${chr}"

    regex="chr${chr}([._/\\-]|$)"

    echo "Copying chr${chr} -> ${DEST_BASE}/chr${chr}"
    grep -E "$regex" "$tmp_pfiles" | \
      gsutil -m -u "$GOOGLE_PROJECT" cp -I "${DEST_BASE}/chr${chr}/"
  done
fi

echo
echo "Done."

# Variant annotation file

VAT="${CDR_STORAGE_PATH}/wgs/short_read/snpindel/aux/vat/vat_complete.bgz.tsv.gz"
OUT="vat_rsids.clean.tsv"

gsutil -u "$GOOGLE_PROJECT" cat "$VAT" \
  | gzip -dc \
  | awk -F'\t' '
      NR==1{
        for(i=1;i<=NF;i++){
          gsub(/^[[:space:]]+|[[:space:]]+$/,"",$i)
          if($i=="contig" || $i=="chrom" || $i=="chr" || $i=="chromosome") c=i
          if($i=="position" || $i=="pos") p=i
          if($i=="dbsnp_rsid" || $i=="rsid" || $i=="dbSNP_rsid") r=i
        }
        if(!c || !p || !r){
          print "no cols with that name" > "/dev/stderr"
          print "Header: " $0 > "/dev/stderr"
          exit 2
        }
        print "contig\tposition\tdbsnp_rsid"
        next
      }
      {
        rs = $r
        if(rs=="" || rs=="." || rs=="NA") next
        if(rs !~ /^rs[0-9]+$/) next
        print $c "\t" $p "\t" rs
      }
    ' > "$OUT"


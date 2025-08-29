#!/bin/bash
set -euo pipefail

###########################################
# Author: Leonard Ndwiga
# Date: 7th August 2025
# Purpose: ONT De Novo Assembly Pipeline (ONT-only)
#          Checks dependencies, gives install commands
###########################################

usage() {
    cat <<'EOF'
ONT De Novo Assembly Pipeline (ONT-only)
=======================================
Description:
  - Concatenate barcode FASTQs
  - Adapter trimming (Porechop)
  - QC (FastQC, MultiQC)
  - Filtering (NanoFilt)
  - Assembly (SPAdes)
  - Evaluation (QUAST)

Supported Reads:
  Only Oxford Nanopore Technologies (ONT) reads are supported.

Usage:
  $(basename "$0") -i <input_dir> -o <output_dir> [-t <threads>]
  $(basename "$0") --input <input_dir> --output <output_dir> [--threads <threads>]

Options:
  -i, --input       Path to input directory containing barcode* folders  (required)
  -o, --output      Path to output directory                              (required)
  -t, --threads     Number of threads to use (default: 4)                 (optional)
  -h, --help        Show help and exit

Example:
  $(basename "$0") -i /data/run1/rawReads -o /work/assembly_out -t 12
EOF
}

INPUT_DIR=""
OUTPUT_DIR=""
THREADS="4"     # Default threads

# ---- Parse arguments ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)   INPUT_DIR="$2"; shift 2;;
    -o|--output)  OUTPUT_DIR="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    -h|--help)    usage; exit 0;;
    *) echo "❌ Unknown option: $1"; usage; exit 1;;
  esac
done

if [[ -z "${INPUT_DIR}" || -z "${OUTPUT_DIR}" ]]; then
  echo "❌ Missing required arguments."; usage; exit 1
fi

# threads must be a positive integer
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
  echo "❌ --threads/-t must be a positive integer. Got: '$THREADS'"; exit 1
fi

# ---- normalize to absolute paths ----
abs_path() {
  local p="$1"
  if [[ -d "$p" ]]; then (cd "$p" && pwd -P)
  else (cd "$(dirname "$p")" && echo "$(pwd -P)/$(basename "$p")")
  fi
}
INPUT_DIR="$(abs_path "$INPUT_DIR")"
OUTPUT_DIR="$(abs_path "$OUTPUT_DIR")"

# ---- Dependency check ----
declare -A REQS=(
  ["porechop"]="porechop"
  ["fastqc"]="fastqc"
  ["multiqc"]="multiqc"
  ["NanoFilt"]="nanofilt"
  ["spades.py"]="spades"
  ["quast.py"]="quast"
)

echo "🔍 Checking required tools..."
missing_pkgs=()

for cmd in "${!REQS[@]}"; do
  if command -v "$cmd" >/dev/null 2>&1; then
    echo "✅ '$cmd' is installed."
  else
    echo "⚠️  Tool '$cmd' is missing."
    missing_pkgs+=("${REQS[$cmd]}")
  fi
done

if (( ${#missing_pkgs[@]} > 0 )); then
  echo ""
  echo "❌ The following required packages are missing:"
  for pkg in "${missing_pkgs[@]}"; do
    echo "    - $pkg"
  done
  echo ""
  echo "Install them with:"
  for pkg in "${missing_pkgs[@]}"; do
    echo "    conda install bioconda::${pkg}"
  done
  exit 1
fi

# ---- create dirs AFTER normalization ----
READS_CAT="${OUTPUT_DIR}/readsCat"
MULTIQC_DIR="${OUTPUT_DIR}/multiqc"
mkdir -p "$READS_CAT" "$MULTIQC_DIR"

# safer globbing
shopt -s nullglob

# timing helpers
SCRIPT_START=$SECONDS
CHUNK_START=$SECONDS
log_time() {
  local label="$1"
  local now=$SECONDS
  local duration=$(( now - CHUNK_START ))
  echo "${label} took ${duration} seconds"
  CHUNK_START=$now
}

# STEP 1: concatenate FASTQs
cd "$INPUT_DIR"
barcode_folders=(barcode*/)

if (( ${#barcode_folders[@]} == 0 )); then
  echo "❌ No barcode folders found in $INPUT_DIR"; exit 1
fi

echo "Found ${#barcode_folders[@]} barcode folders:"
printf ' - %s\n' "${barcode_folders[@]}"

for folder in "${barcode_folders[@]}"; do
  folder_name="${folder%/}"
  output_file="${READS_CAT}/${folder_name}.fastq.gz"

  echo "Combining chunks in $folder_name..."
  fq_files=("${folder}"/*.fastq.gz)
  if (( ${#fq_files[@]} == 0 )); then
    echo "  ⚠️ No FASTQ files in $folder_name. Skipping."
    continue
  fi

  mkdir -p "$(dirname "$output_file")"
  zcat "${fq_files[@]}" | gzip > "$output_file"
  echo "  → Created: ${output_file}"
done
log_time "Concatenate barcode FASTQs (zcat)"

# STEP 2: Porechop + FastQC
cd "$READS_CAT"
barcode_fastqs=(*.fastq.gz)
if (( ${#barcode_fastqs[@]} == 0 )); then
  echo "❌ No .fastq.gz files found in $READS_CAT"; exit 1
fi

for fq in "${barcode_fastqs[@]}"; do
  sample="${fq%.fastq.gz}"
  echo "Adapter removal + QC for $sample ..."
  sample_dir="${OUTPUT_DIR}/${sample}"
  mkdir -p "$sample_dir/fastqc"

  porechop --threads "$THREADS" -i "$fq" -o "$sample_dir/${sample}_porechop.fastq.gz"
  log_time "$sample - Porechop"

  fastqc --threads "$THREADS" -o "$sample_dir/fastqc" "$sample_dir/${sample}_porechop.fastq.gz"
  log_time "$sample - FastQC"
done

# STEP 3: MultiQC
echo "Running MultiQC..."
multiqc "$OUTPUT_DIR" -o "$MULTIQC_DIR"
log_time "MultiQC summary"

# STEP 4–6: NanoFilt → SPAdes → QUAST
for fq in "${barcode_fastqs[@]}"; do
  sample="${fq%.fastq.gz}"
  echo "🔧 Continuing with $sample (Trimming → Assembly → Evaluation)..."
  sample_dir="${OUTPUT_DIR}/${sample}"

  zcat "$sample_dir/${sample}_porechop.fastq.gz" | NanoFilt -q 10 -l 500 | gzip > "$sample_dir/${sample}_trimmed.fastq.gz"
  log_time "$sample - Trimming (NanoFilt)"

  spades.py --threads "$THREADS" --isolate -s "$sample_dir/${sample}_trimmed.fastq.gz" -o "$sample_dir/spades_output"
  log_time "$sample - SPAdes Assembly"

  quast.py -t "$THREADS" "$sample_dir/spades_output/contigs.fasta" -o "$sample_dir/quast_output"
  log_time "$sample - QUAST Evaluation"

  echo "✅ $sample complete. Contigs: $sample_dir/spades_output/contigs.fasta"
done

TOTAL_RUNTIME=$((SECONDS - SCRIPT_START))
echo "Pipeline complete."
echo "Total script runtime: ${TOTAL_RUNTIME} seconds"

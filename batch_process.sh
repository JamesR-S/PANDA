#!/usr/bin/env bash

set -euo pipefail

# Defaults

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT_LIST=""
OUTPUT_PREFIX=""

THREADS=8
CLASSIFICATION_THRESHOLD=0.55

RESOURCE_DIR="${SCRIPT_DIR}/resources"
BIN_DIR="${SCRIPT_DIR}/bin"

KEEP_TEMP=false
SETUP_ONLY=false
FORCE_DOWNLOAD=false

# Replace this once your Zenodo record is final.
CUSTOM_RESOURCES_URL="https://zenodo.org/record/21217549/files/resources.tar.gz"

# Replace these if the official URLs differ.
LASER_URL="https://csg.sph.umich.edu/chaolong/LASER/LASER-2.04.tar.gz"
FASTA_URL="ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"

REFERENCE_FASTA="${RESOURCE_DIR}/reference/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
REFERENCE_FASTA_INDEX="${REFERENCE_FASTA}.fai"
REFERENCE_BED="${RESOURCE_DIR}/reference/New_Ref_800k.bed"
REFERENCE_SITE="${RESOURCE_DIR}/reference/New_Ref_800k.site"
REFERENCE_GENO="${RESOURCE_DIR}/reference/New_Ref_800k.geno"
REFERENCE_COORD="${RESOURCE_DIR}/reference/New_Ref_800k.RefPC.coord"

LASER_DIR="${BIN_DIR}/LASER-2.04"
LASER_EXE="${LASER_DIR}/laser"
PILEUP2SEQ="${LASER_DIR}/pileup2seq/pileup2seq.py"
CLASSIFICATION_SCRIPT="${BIN_DIR}/classification.rscript"

TEMP_DIR=""
PILEUP_DIR=""

# Helper functions

info() {
    echo "[INFO] $*"
}

warn() {
    echo "[WARNING] $*" >&2
}

error() {
    echo "[ERROR] $*" >&2
}

die() {
    error "$*"
    exit 1
}

usage() {
    cat <<EOF
Usage:
  $(basename "$0") --input bam_list.txt --output output_prefix [options]

Required:
  -i, --input FILE              Text file containing one BAM/CRAM path per line
  -o, --output PREFIX           Output prefix

Options:
  -t, --threads INT             Number of parallel jobs to use [default: ${THREADS}]
      --threshold FLOAT         Classification probability threshold [default: ${CLASSIFICATION_THRESHOLD}]

      --resource-dir DIR        Resource directory [default: ${RESOURCE_DIR}]
      --bin-dir DIR             Binary/software directory [default: ${BIN_DIR}]

      --custom-resources-url URL
                                URL for custom resources tar.gz archive
      --laser-url URL           URL for LASER tar.gz archive
      --fasta-url URL           URL for reference FASTA .gz file

      --setup-only              Download/check resources and exit
      --force-download          Re-download resources even if files already exist
      --keep-temp               Keep temporary files after successful run

  -h, --help                    Show this help message

Examples:
  $(basename "$0") \\
    --input bams.txt \\
    --output cohort1 \\
    --threads 16

  $(basename "$0") --setup-only

Outputs:
  PREFIX.coord
  PREFIX.population.classification.txt
  unprocessed_ids.txt

EOF
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

download_file() {
    local url="$1"
    local output="$2"

    if [[ -f "$output" && "$FORCE_DOWNLOAD" == false ]]; then
        info "File already exists, skipping download: $output"
        return
    fi

    mkdir -p "$(dirname "$output")"

    info "Downloading:"
    info "  $url"
    info "to:"
    info "  $output"

    if command_exists curl; then
        curl -L --fail -o "$output" "$url"
    elif command_exists wget; then
        wget -O "$output" "$url"
    else
        die "Neither curl nor wget is installed. Please install one of them, or download files manually."
    fi
}

check_required_commands() {
    local missing=()

    for cmd in samtools parallel Rscript python2 tar gzip find split head tail wc awk; do
        if ! command_exists "$cmd"; then
            missing+=("$cmd")
        fi
    done

    if [[ "${#missing[@]}" -gt 0 ]]; then
        error "Missing required command-line tool(s): ${missing[*]}"
        echo
        echo "Please install these before running the pipeline."
        echo
        echo "Typical installation examples:"
        echo "  conda install -c bioconda samtools parallel"
        echo "  conda install -c conda-forge r-base"
        echo "  conda install -c conda-forge python=2.7"
        echo
        exit 1
    fi
}

check_required_files() {
    local missing=()

    for file in \
        "$REFERENCE_FASTA" \
        "$REFERENCE_BED" \
        "$REFERENCE_SITE" \
        "$REFERENCE_GENO" \
        "$REFERENCE_COORD" \
        "$LASER_EXE" \
        "$PILEUP2SEQ" \
        "$CLASSIFICATION_SCRIPT"
    do
        if [[ ! -e "$file" ]]; then
            missing+=("$file")
        fi
    done

    if [[ "${#missing[@]}" -gt 0 ]]; then
        error "Missing required file(s):"
        printf '  %s\n' "${missing[@]}"
        echo
        echo "You can try to download/setup resources using:"
        echo
        echo "  $(basename "$0") --setup-only"
        echo
        echo "Or provide alternative locations with:"
        echo
        echo "  --resource-dir /path/to/resources --bin-dir /path/to/bin"
        echo
        exit 1
    fi
}

download_and_prepare_resources() {
    mkdir -p "$RESOURCE_DIR" "$BIN_DIR"

    local custom_archive="${SCRIPT_DIR}/resources.tar.gz"
    local laser_archive="${SCRIPT_DIR}/LASER-2.04.tar.gz"
    local fasta_gz="${REFERENCE_FASTA}.gz"

    info "Checking custom resources."

    if [[ "$FORCE_DOWNLOAD" == true || ! -f "$REFERENCE_BED" || ! -f "$REFERENCE_SITE" || ! -f "$REFERENCE_GENO" || ! -f "$REFERENCE_COORD" ]]; then
        download_file "$CUSTOM_RESOURCES_URL" "$custom_archive"

        info "Extracting custom resources archive."

        mkdir -p "$RESOURCE_DIR"

        # This assumes the archive contains either:
        #   resources/reference/...
        # or files that extract into the expected resources directory.
        tar -xzf "$custom_archive" -C "$SCRIPT_DIR"

        if [[ ! -d "$RESOURCE_DIR" ]]; then
            die "Custom resources archive was extracted, but ${RESOURCE_DIR} was not found."
        fi
    else
        info "Custom resources appear to be present."
    fi

    info "Checking LASER."

    if [[ "$FORCE_DOWNLOAD" == true || ! -x "$LASER_EXE" || ! -f "$PILEUP2SEQ" ]]; then
        download_file "$LASER_URL" "$laser_archive"

        info "Extracting LASER archive."

        mkdir -p "$BIN_DIR"
        tar -xzf "$laser_archive" -C "$BIN_DIR"

        if [[ -f "$LASER_EXE" && ! -x "$LASER_EXE" ]]; then
            chmod +x "$LASER_EXE" 2>/dev/null || true
        fi
    else
        info "LASER appears to be present."
    fi

    info "Checking reference FASTA."

    if [[ "$FORCE_DOWNLOAD" == true || ! -f "$REFERENCE_FASTA" ]]; then
        mkdir -p "$(dirname "$REFERENCE_FASTA")"

        download_file "$FASTA_URL" "$fasta_gz"

        info "Decompressing FASTA."
        gzip -f -d "$fasta_gz"
    else
        info "Reference FASTA appears to be present."
    fi
}

prepare_reference_index() {
    if [[ -f "$REFERENCE_FASTA_INDEX" ]]; then
        info "Reference FASTA index found: $REFERENCE_FASTA_INDEX"
    else
        info "Reference FASTA index not found. Creating it with samtools faidx."
        samtools faidx "$REFERENCE_FASTA"
    fi

    if [[ ! -f "$REFERENCE_FASTA_INDEX" ]]; then
        die "Failed to create FASTA index: $REFERENCE_FASTA_INDEX"
    fi
}

cleanup() {
    if [[ -n "${TEMP_DIR:-}" && -d "$TEMP_DIR" && "$KEEP_TEMP" == false ]]; then
        info "Cleaning up temporary directory: $TEMP_DIR"
        rm -rf "$TEMP_DIR"
    elif [[ -n "${TEMP_DIR:-}" && -d "$TEMP_DIR" ]]; then
        info "Keeping temporary directory: $TEMP_DIR"
    fi
}

validate_numeric_options() {
    [[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be a positive integer."
    [[ "$THREADS" -gt 0 ]] || die "--threads must be greater than zero."

    awk "BEGIN { exit !($CLASSIFICATION_THRESHOLD >= 0 && $CLASSIFICATION_THRESHOLD <= 1) }" \
        || die "--threshold must be a number between 0 and 1."
}

# Parse arguments

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT_LIST="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        --threshold)
            CLASSIFICATION_THRESHOLD="$2"
            shift 2
            ;;
        --resource-dir)
            RESOURCE_DIR="$2"
            shift 2
            ;;
        --bin-dir)
            BIN_DIR="$2"
            shift 2
            ;;
        --custom-resources-url)
            CUSTOM_RESOURCES_URL="$2"
            shift 2
            ;;
        --laser-url)
            LASER_URL="$2"
            shift 2
            ;;
        --fasta-url)
            FASTA_URL="$2"
            shift 2
            ;;
        --setup-only)
            SETUP_ONLY=true
            shift
            ;;
        --force-download)
            FORCE_DOWNLOAD=true
            shift
            ;;
        --keep-temp)
            KEEP_TEMP=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "Unknown option: $1. Use --help for usage."
            ;;
    esac
done

# Recompute paths after option parsing

REFERENCE_FASTA="${RESOURCE_DIR}/reference/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
REFERENCE_FASTA_INDEX="${REFERENCE_FASTA}.fai"
REFERENCE_BED="${RESOURCE_DIR}/reference/New_Ref_800k.bed"
REFERENCE_SITE="${RESOURCE_DIR}/reference/New_Ref_800k.site"
REFERENCE_GENO="${RESOURCE_DIR}/reference/New_Ref_800k.geno"
REFERENCE_COORD="${RESOURCE_DIR}/reference/New_Ref_800k.RefPC.coord"

LASER_DIR="${BIN_DIR}/LASER-2.04"
LASER_EXE="${LASER_DIR}/laser"
PILEUP2SEQ="${LASER_DIR}/pileup2seq/pileup2seq.py"
CLASSIFICATION_SCRIPT="${BIN_DIR}/classification.rscript"

validate_numeric_options
check_required_commands

download_and_prepare_resources

if [[ "$SETUP_ONLY" == true ]]; then
    check_required_files
    prepare_reference_index
    info "Setup complete."
    exit 0
fi

# Validate required run inputs

[[ -n "$INPUT_LIST" ]] || die "Missing required argument: --input"
[[ -n "$OUTPUT_PREFIX" ]] || die "Missing required argument: --output"

[[ -f "$INPUT_LIST" ]] || die "Input list does not exist: $INPUT_LIST"
[[ -s "$INPUT_LIST" ]] || die "Input list is empty: $INPUT_LIST"

check_required_files
prepare_reference_index

# Pipeline setup

TEMP_DIR="${SCRIPT_DIR}/temp_files_${OUTPUT_PREFIX}_$$"
PILEUP_DIR="${TEMP_DIR}/pileup_files"

mkdir -p "$PILEUP_DIR"

trap cleanup EXIT

export SCRIPT_DIR
export RESOURCE_DIR
export BIN_DIR
export REFERENCE_FASTA
export REFERENCE_BED
export REFERENCE_SITE
export REFERENCE_GENO
export REFERENCE_COORD
export LASER_EXE
export PILEUP2SEQ
export PILEUP_DIR
export TEMP_DIR

info "Pipeline configuration:"
info "  Script directory:              $SCRIPT_DIR"
info "  Input list:                    $INPUT_LIST"
info "  Output prefix:                 $OUTPUT_PREFIX"
info "  Resource directory:            $RESOURCE_DIR"
info "  Binary directory:              $BIN_DIR"
info "  Parallel jobs:                 $THREADS"
info "  LASER mode:                    single-threaded per .seq file"
info "  Classification threshold:      $CLASSIFICATION_THRESHOLD"
info "  Temporary directory:           $TEMP_DIR"

# Step 1: Make pileup files

info "Step 1/5: Making pileup files."

parallel \
    --halt soon,fail=1 \
    -a "$INPUT_LIST" \
    -j "$THREADS" \
    --joblog "${OUTPUT_PREFIX}.make_pileup.log" \
    'samtools mpileup -q 30 -Q 20 -f "$REFERENCE_FASTA" -l "$REFERENCE_BED" "{}" > "$PILEUP_DIR/{/.}.pileup"'

info "Pileup files created."

# Step 2: Convert pileup files to .seq

info "Step 2/5: Creating .seq files from pileup files."

PILEUP_COUNT=$(find "$PILEUP_DIR" -type f -name "*.pileup" | wc -l | tr -d ' ')

if [[ "$PILEUP_COUNT" -eq 0 ]]; then
    die "No pileup files were created."
fi

# Create approximately one pileup2seq batch per thread.
# Example:
#   10,000 pileup files / 20 threads = 500 files per batch
BATCH_SIZE=$(( (PILEUP_COUNT + THREADS - 1) / THREADS ))

info "Found $PILEUP_COUNT pileup file(s)."
info "Creating up to $THREADS pileup2seq batch(es), with up to $BATCH_SIZE pileup file(s) per batch."

find "$PILEUP_DIR" -type f -name "*.pileup" | sort | split -l "$BATCH_SIZE" - "${TEMP_DIR}/temp_batch_"

parallel \
    --halt soon,fail=1 \
    -j "$THREADS" \
    'batch=$(cat "{}"); python2 "$PILEUP2SEQ" -m "$REFERENCE_SITE" -o "{.}_800k" $batch; rm "{}"' \
    ::: "${TEMP_DIR}"/temp_batch_*

info ".seq files created."

# Step 3: Run LASER

info "Step 3/5: Running LASER on each .seq file."

find "$TEMP_DIR" -name "*.seq" | sort | parallel \
    --halt soon,fail=1 \
    -j "$THREADS" \
    '"$LASER_EXE" -nt 1 -s "{}" -g "$REFERENCE_GENO" -c "$REFERENCE_COORD" -o "{.}" -k 10'

info ".coord files created."

# Step 4: Concatenate .coord files

info "Step 4/5: Concatenating .coord files."

mapfile -t coord_files < <(find "$TEMP_DIR" -name "*.coord" | sort)

if [[ "${#coord_files[@]}" -eq 0 ]]; then
    die "No .coord files were produced by LASER."
fi

head -n 1 "${coord_files[0]}" > "${OUTPUT_PREFIX}.coord"

for coord_file in "${coord_files[@]}"; do
    tail -n +2 "$coord_file" >> "${OUTPUT_PREFIX}.coord"
done

info "Concatenated coordinates written to: ${OUTPUT_PREFIX}.coord"

# Step 5: Classification

info "Step 5/5: Running classification model."

Rscript "$CLASSIFICATION_SCRIPT" \
    "${OUTPUT_PREFIX}.coord" \
    "$CLASSIFICATION_THRESHOLD" \
    "${OUTPUT_PREFIX}.population.classification.txt"

info "Classification written to: ${OUTPUT_PREFIX}.population.classification.txt"
info "Processing complete."

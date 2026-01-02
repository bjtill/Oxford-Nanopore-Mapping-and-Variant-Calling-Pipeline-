#!/bin/bash

################################################################################
# Script: 04a_call_variants_clair3.sh
# Description: Call variants using Clair3 (deep learning-based caller)
# Usage: ./04a_call_variants_clair3.sh <input_bam> <reference> <output_dir> <model> <threads>
################################################################################

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <input_bam> <reference_fasta> <output_dir> <model> <threads>"
    echo ""
    echo "Arguments:"
    echo "  input_bam        Path to sorted and indexed BAM file"
    echo "  reference_fasta  Path to reference genome FASTA file"
    echo "  output_dir       Directory for Clair3 output"
    echo "  model            Clair3 model name:"
    echo "                   - r1041_e82_400bps_sup_v500 (R10.4.1, recommended)"
    echo "                   - r1041_e82_400bps_sup_v410 (R10.4.1, older)"
    echo "                   - r941_prom_sup_g5014 (R9.4.1)"
    echo "                   - r941_prom_hac_g360+g422 (R9.4.1)"
    echo "  threads          Number of threads to use"
    echo ""
    echo "Example:"
    echo "  $0 aligned.sorted.bam reference.fasta clair3_output r1041_e82_400bps_sup_v500 16"
    exit 1
}

# Check arguments
if [ $# -ne 5 ]; then
    usage
fi

INPUT_BAM="$1"
REFERENCE="$2"
OUTPUT_DIR="$3"
MODEL="$4"
THREADS="$5"

# Validate inputs
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file not found: $INPUT_BAM"
    exit 1
fi

if [ ! -f "${INPUT_BAM}.bai" ]; then
    echo "Error: BAM index not found. Please run 03_process_bam.sh first."
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference FASTA file not found: $REFERENCE"
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error: Threads must be a positive integer"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_DIR}/clair3.log"

echo "========================================" | tee "$LOG_FILE"
echo "Variant Calling with Clair3" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "Input BAM: $INPUT_BAM" | tee -a "$LOG_FILE"
echo "Reference: $REFERENCE" | tee -a "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Model: $MODEL" | tee -a "$LOG_FILE"
echo "Threads: $THREADS" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Activate Clair3 conda environment
echo "Activating conda environment: clair3" | tee -a "$LOG_FILE"
eval "$(conda shell.bash hook)"
conda activate clair3

# Check if Clair3 is available
if ! command -v run_clair3.sh &> /dev/null; then
    echo "Error: Clair3 not found. Please install it in the clair3 conda environment."
    exit 1
fi

# Determine models directory based on where run_clair3.sh is located
CLAIR3_BIN=$(dirname $(which run_clair3.sh))
MODEL_BASE="${CLAIR3_BIN}/models"

# Construct model path
MODEL_PATH="${MODEL_BASE}/${MODEL}"

# Check if model exists
if [ ! -d "$MODEL_PATH" ]; then
    echo "Error: Model not found at: $MODEL_PATH" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    echo "Available models in ${MODEL_BASE}:" | tee -a "$LOG_FILE"
    ls -1 "${MODEL_BASE}/" 2>/dev/null | grep -v "^ont" | tee -a "$LOG_FILE" || echo "No models directory found"
    echo "" | tee -a "$LOG_FILE"
    echo "Please use one of the available models."
    exit 1
fi

echo "Using model path: $MODEL_PATH" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Run Clair3
echo "Running Clair3 variant calling..." | tee -a "$LOG_FILE"
echo "This may take several hours depending on genome size and coverage..." | tee -a "$LOG_FILE"

run_clair3.sh \
    --bam_fn="$INPUT_BAM" \
    --ref_fn="$REFERENCE" \
    --output="$OUTPUT_DIR" \
    --threads="$THREADS" \
    --platform=ont \
    --model_path="$MODEL_PATH" \
    --include_all_ctgs \
    --no_phasing_for_fa 2>&1 | tee -a "$LOG_FILE"

if [ $? -eq 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "Clair3 variant calling completed successfully!" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    echo "Output files:" | tee -a "$LOG_FILE"
    echo "  - ${OUTPUT_DIR}/merge_output.vcf.gz (final VCF)" | tee -a "$LOG_FILE"
    echo "  - ${OUTPUT_DIR}/merge_output.vcf.gz.tbi (VCF index)" | tee -a "$LOG_FILE"
    
    # Quick variant count
    if [ -f "${OUTPUT_DIR}/merge_output.vcf.gz" ]; then
        VARIANT_COUNT=$(zgrep -v "^#" "${OUTPUT_DIR}/merge_output.vcf.gz" | wc -l)
        echo "" | tee -a "$LOG_FILE"
        echo "Total variants called: $VARIANT_COUNT" | tee -a "$LOG_FILE"
    fi
else
    echo "Error: Clair3 variant calling failed" | tee -a "$LOG_FILE"
    exit 1
fi

echo "" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

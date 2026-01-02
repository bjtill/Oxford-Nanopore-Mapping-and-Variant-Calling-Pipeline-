#!/bin/bash

################################################################################
# Script: 04b_call_variants_pepper.sh
# Description: Call variants using PEPPER-Margin-DeepVariant
# Usage: ./04b_call_variants_pepper.sh <input_bam> <reference> <output_dir> <model> <threads>
################################################################################

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <input_bam> <reference_fasta> <output_dir> <model> <threads>"
    echo ""
    echo "Arguments:"
    echo "  input_bam        Path to sorted and indexed BAM file"
    echo "  reference_fasta  Path to reference genome FASTA file"
    echo "  output_dir       Directory for PEPPER output"
    echo "  model            Sequencing model:"
    echo "                   - ont_r9_guppy5_sup"
    echo "                   - ont_r10_guppy5_sup"
    echo "                   - ont_r10_dorado_5khz_sup (recommended for recent R10 data)"
    echo "  threads          Number of threads to use"
    echo ""
    echo "Example:"
    echo "  $0 aligned.sorted.bam reference.fasta pepper_output ont_r10_dorado_5khz_sup 16"
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

# Validate model
case "$MODEL" in
    ont_r9_guppy5_sup|ont_r10_guppy5_sup|ont_r10_dorado_5khz_sup)
        ;;
    *)
        echo "Error: Invalid model: $MODEL"
        echo "Valid options: ont_r9_guppy5_sup, ont_r10_guppy5_sup, ont_r10_dorado_5khz_sup"
        exit 1
        ;;
esac

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_DIR}/pepper.log"

echo "========================================" | tee "$LOG_FILE"
echo "Variant Calling with PEPPER-Margin-DeepVariant" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "Input BAM: $INPUT_BAM" | tee -a "$LOG_FILE"
echo "Reference: $REFERENCE" | tee -a "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Model: $MODEL" | tee -a "$LOG_FILE"
echo "Threads: $THREADS" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Activate PEPPER conda environment
echo "Activating conda environment: pepper" | tee -a "$LOG_FILE"
eval "$(conda shell.bash hook)"
conda activate pepper

# Check if PEPPER is available
if ! command -v run_pepper_margin_deepvariant &> /dev/null; then
    echo "Error: PEPPER-Margin-DeepVariant not found."
    echo "Please install it in the pepper conda environment."
    echo "Installation: pip install pepper-deepvariant-runner"
    exit 1
fi

# Run PEPPER-Margin-DeepVariant
echo "Running PEPPER-Margin-DeepVariant variant calling..." | tee -a "$LOG_FILE"
echo "This may take several hours depending on genome size and coverage..." | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

run_pepper_margin_deepvariant call_variant \
    -b "$INPUT_BAM" \
    -f "$REFERENCE" \
    -o "$OUTPUT_DIR" \
    -t "$THREADS" \
    --${MODEL} \
    --phased_output 2>&1 | tee -a "$LOG_FILE"

if [ $? -eq 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "PEPPER-Margin-DeepVariant variant calling completed successfully!" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    echo "Output files:" | tee -a "$LOG_FILE"
    echo "  - ${OUTPUT_DIR}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz" | tee -a "$LOG_FILE"
    echo "  - ${OUTPUT_DIR}/PEPPER_MARGIN_DEEPVARIANT_PHASED_FINAL_OUTPUT.vcf.gz" | tee -a "$LOG_FILE"
    
    # Quick variant count
    if [ -f "${OUTPUT_DIR}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz" ]; then
        VARIANT_COUNT=$(zgrep -v "^#" "${OUTPUT_DIR}/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz" | wc -l)
        echo "" | tee -a "$LOG_FILE"
        echo "Total variants called: $VARIANT_COUNT" | tee -a "$LOG_FILE"
    fi
else
    echo "Error: PEPPER-Margin-DeepVariant variant calling failed" | tee -a "$LOG_FILE"
    exit 1
fi

echo "" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

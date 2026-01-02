#!/bin/bash

################################################################################
# Script: 01_qc.sh
# Description: Quality control of Oxford Nanopore reads using NanoPlot
# Usage: ./01_qc.sh <input_fastq> <output_dir> <threads>
################################################################################

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <input_fastq> <output_dir> <threads>"
    echo ""
    echo "Arguments:"
    echo "  input_fastq   Path to input FASTQ file (can be gzipped)"
    echo "  output_dir    Directory for QC output"
    echo "  threads       Number of threads to use"
    echo ""
    echo "Example:"
    echo "  $0 reads.fastq.gz qc_results 16"
    exit 1
}

# Check arguments
if [ $# -ne 3 ]; then
    usage
fi

INPUT_FASTQ="$1"
OUTPUT_DIR="$2"
THREADS="$3"

# Validate inputs
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "Error: Input FASTQ file not found: $INPUT_FASTQ"
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error: Threads must be a positive integer"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_DIR}/01_qc.log"

echo "========================================" | tee "$LOG_FILE"
echo "Oxford Nanopore QC with NanoPlot" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "Input FASTQ: $INPUT_FASTQ" | tee -a "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Threads: $THREADS" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Activate conda environment
echo "Activating conda environment: nanopore_pipeline" | tee -a "$LOG_FILE"
eval "$(conda shell.bash hook)"
conda activate nanopore_pipeline

# Check if NanoPlot is available
if ! command -v NanoPlot &> /dev/null; then
    echo "Error: NanoPlot not found. Please install it in the nanopore_pipeline conda environment."
    exit 1
fi

# Run NanoPlot
echo "Running NanoPlot..." | tee -a "$LOG_FILE"
NanoPlot \
    --fastq "$INPUT_FASTQ" \
    --outdir "$OUTPUT_DIR" \
    --threads "$THREADS" \
    --plots dot \
    --legacy hex \
    --N50 \
    --loglength 2>&1 | tee -a "$LOG_FILE"

if [ $? -eq 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "QC completed successfully!" | tee -a "$LOG_FILE"
    echo "Results saved to: $OUTPUT_DIR" | tee -a "$LOG_FILE"
    echo "View NanoPlot-report.html for detailed results" | tee -a "$LOG_FILE"
else
    echo "Error: NanoPlot failed" | tee -a "$LOG_FILE"
    exit 1
fi

echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

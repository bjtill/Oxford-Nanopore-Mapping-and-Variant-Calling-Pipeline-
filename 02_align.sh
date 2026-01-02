#!/bin/bash

################################################################################
# Script: 02_align.sh
# Description: Align Oxford Nanopore reads to reference genome using minimap2
# Usage: ./02_align.sh <input_fastq> <reference_fasta> <output_bam> <threads>
################################################################################

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <input_fastq> <reference_fasta> <output_bam> <threads>"
    echo ""
    echo "Arguments:"
    echo "  input_fastq      Path to input FASTQ file (can be gzipped)"
    echo "  reference_fasta  Path to reference genome FASTA file"
    echo "  output_bam       Path for output BAM file"
    echo "  threads          Number of threads to use"
    echo ""
    echo "Example:"
    echo "  $0 reads.fastq.gz reference.fasta aligned.bam 16"
    exit 1
}

# Check arguments
if [ $# -ne 4 ]; then
    usage
fi

INPUT_FASTQ="$1"
REFERENCE="$2"
OUTPUT_BAM="$3"
THREADS="$4"

# Validate inputs
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "Error: Input FASTQ file not found: $INPUT_FASTQ"
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

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_BAM")
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_BAM%.bam}.align.log"

echo "========================================" | tee "$LOG_FILE"
echo "Oxford Nanopore Read Alignment" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "Input FASTQ: $INPUT_FASTQ" | tee -a "$LOG_FILE"
echo "Reference: $REFERENCE" | tee -a "$LOG_FILE"
echo "Output BAM: $OUTPUT_BAM" | tee -a "$LOG_FILE"
echo "Threads: $THREADS" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Activate conda environment
echo "Activating conda environment: nanopore_pipeline" | tee -a "$LOG_FILE"
eval "$(conda shell.bash hook)"
conda activate nanopore_pipeline

# Check if tools are available
if ! command -v minimap2 &> /dev/null; then
    echo "Error: minimap2 not found"
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found"
    exit 1
fi

# Index reference if not already done
if [ ! -f "${REFERENCE}.fai" ]; then
    echo "Indexing reference genome..." | tee -a "$LOG_FILE"
    samtools faidx "$REFERENCE" 2>&1 | tee -a "$LOG_FILE"
fi

# Run minimap2 alignment
# Using -ax map-ont for Oxford Nanopore reads
# -a: output SAM format
# -x map-ont: preset for Oxford Nanopore vs reference mapping
# -t: number of threads
# -y: copy read tags to output
# --MD: output MD tag for variant calling
# --eqx: use =/X CIGAR operators for match/mismatch
echo "Aligning reads with minimap2..." | tee -a "$LOG_FILE"
minimap2 \
    -ax map-ont \
    -t "$THREADS" \
    -y \
    --MD \
    --eqx \
    "$REFERENCE" \
    "$INPUT_FASTQ" 2>> "$LOG_FILE" | \
samtools view \
    -@ 4 \
    -b \
    -h \
    -o "$OUTPUT_BAM" \
    - 2>> "$LOG_FILE"

if [ $? -eq 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "Alignment completed successfully!" | tee -a "$LOG_FILE"
    echo "Output BAM: $OUTPUT_BAM" | tee -a "$LOG_FILE"
    
    # Quick alignment statistics
    echo "" | tee -a "$LOG_FILE"
    echo "Quick alignment statistics:" | tee -a "$LOG_FILE"
    samtools flagstat "$OUTPUT_BAM" 2>&1 | tee -a "$LOG_FILE"
else
    echo "Error: Alignment failed" | tee -a "$LOG_FILE"
    exit 1
fi

echo "" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

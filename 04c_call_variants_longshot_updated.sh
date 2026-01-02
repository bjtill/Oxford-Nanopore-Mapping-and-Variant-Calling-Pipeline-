#!/bin/bash

################################################################################
# Script: 04c_call_variants_longshot.sh
# Description: Call variants using Longshot (fast SNV caller for long reads)
# Usage: ./04c_call_variants_longshot.sh <input_bam> <reference> <output_vcf> <threads>
################################################################################

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <input_bam> <reference_fasta> <output_vcf> <threads>"
    echo ""
    echo "Arguments:"
    echo "  input_bam        Path to sorted and indexed BAM file"
    echo "  reference_fasta  Path to reference genome FASTA file"
    echo "  output_vcf       Path for output VCF file"
    echo "  threads          Number of threads to use"
    echo ""
    echo "Example:"
    echo "  $0 aligned.sorted.bam reference.fasta longshot_variants.vcf 16"
    echo ""
    echo "Note: Longshot is optimized for diploid SNV calling."
    echo "      For indels, consider using Clair3."
    exit 1
}

# Check arguments
if [ $# -ne 4 ]; then
    usage
fi

INPUT_BAM="$1"
REFERENCE="$2"
OUTPUT_VCF="$3"
THREADS="$4"

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

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_VCF%.vcf}.log"

echo "========================================" | tee "$LOG_FILE"
echo "Variant Calling with Longshot" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "Input BAM: $INPUT_BAM" | tee -a "$LOG_FILE"
echo "Reference: $REFERENCE" | tee -a "$LOG_FILE"
echo "Output VCF: $OUTPUT_VCF" | tee -a "$LOG_FILE"
echo "Threads: $THREADS" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Activate conda environment
echo "Activating conda environment: nanopore_pipeline" | tee -a "$LOG_FILE"
eval "$(conda shell.bash hook)"
conda activate nanopore_pipeline

# Check if Longshot is available
if ! command -v longshot &> /dev/null; then
    echo "Error: Longshot not found. Please install it in the nanopore_pipeline conda environment."
    exit 1
fi

# Check Longshot version to determine available options
echo "Checking Longshot version..." | tee -a "$LOG_FILE"
longshot --version 2>&1 | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Index reference if needed
if [ ! -f "${REFERENCE}.fai" ]; then
    echo "Indexing reference genome..." | tee -a "$LOG_FILE"
    samtools faidx "$REFERENCE" 2>&1 | tee -a "$LOG_FILE"
fi

# Run Longshot with basic parameters that work across versions
echo "Running Longshot variant calling..." | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Basic Longshot command compatible with most versions
longshot \
    --bam "$INPUT_BAM" \
    --ref "$REFERENCE" \
    --out "$OUTPUT_VCF" \
    --min_cov 5 \
    --min_mapq 30 \
    --min_allele_qual 7.0 2>&1 | tee -a "$LOG_FILE"

LONGSHOT_EXIT=$?

if [ $LONGSHOT_EXIT -eq 0 ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "Longshot variant calling completed successfully!" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    
    # Check if VCF was created
    if [ ! -f "$OUTPUT_VCF" ]; then
        echo "Warning: Output VCF not found, but Longshot reported success" | tee -a "$LOG_FILE"
        exit 1
    fi
    
    # Compress and index VCF
    echo "Compressing and indexing VCF..." | tee -a "$LOG_FILE"
    
    # Check if output is already compressed
    if [[ "$OUTPUT_VCF" != *.gz ]]; then
        bgzip -f "$OUTPUT_VCF"
        OUTPUT_VCF="${OUTPUT_VCF}.gz"
    fi
    
    tabix -f -p vcf "$OUTPUT_VCF"
    
    echo "Output files:" | tee -a "$LOG_FILE"
    echo "  - ${OUTPUT_VCF} (compressed VCF)" | tee -a "$LOG_FILE"
    echo "  - ${OUTPUT_VCF}.tbi (VCF index)" | tee -a "$LOG_FILE"
    
    # Quick variant count
    if [ -f "$OUTPUT_VCF" ]; then
        VARIANT_COUNT=$(zgrep -v "^#" "$OUTPUT_VCF" | wc -l)
        SNP_COUNT=$(zgrep -v "^#" "$OUTPUT_VCF" | awk 'length($4)==1 && length($5)==1' | wc -l)
        
        echo "" | tee -a "$LOG_FILE"
        echo "Variant counts:" | tee -a "$LOG_FILE"
        echo "  Total variants: $VARIANT_COUNT" | tee -a "$LOG_FILE"
        echo "  SNPs (estimated): $SNP_COUNT" | tee -a "$LOG_FILE"
    fi
else
    echo "Error: Longshot variant calling failed with exit code $LONGSHOT_EXIT" | tee -a "$LOG_FILE"
    exit 1
fi

echo "" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

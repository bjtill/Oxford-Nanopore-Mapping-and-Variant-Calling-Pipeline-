#!/bin/bash

################################################################################
# Script: 05_filter_variants.sh
# Description: Filter variants based on quality and depth
# Usage: ./05_filter_variants.sh <input_vcf> <output_vcf> <min_qual> <min_depth>
################################################################################

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <input_vcf> <output_vcf> <min_qual> <min_depth>"
    echo ""
    echo "Arguments:"
    echo "  input_vcf    Path to input VCF file (can be gzipped)"
    echo "  output_vcf   Path for filtered output VCF file"
    echo "  min_qual     Minimum variant quality score (QUAL)"
    echo "  min_depth    Minimum read depth (DP)"
    echo ""
    echo "Example:"
    echo "  $0 variants.vcf.gz filtered.vcf.gz 20 5"
    echo ""
    echo "Recommended filtering thresholds:"
    echo "  - QUAL >= 20 (Phred-scaled quality)"
    echo "  - DP >= 5 (minimum coverage)"
    echo "  - For high-confidence calls: QUAL >= 30, DP >= 10"
    exit 1
}

# Check arguments
if [ $# -ne 4 ]; then
    usage
fi

INPUT_VCF="$1"
OUTPUT_VCF="$2"
MIN_QUAL="$3"
MIN_DEPTH="$4"

# Validate inputs
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF file not found: $INPUT_VCF"
    exit 1
fi

if ! [[ "$MIN_QUAL" =~ ^[0-9]+\.?[0-9]*$ ]]; then
    echo "Error: Minimum quality must be a number"
    exit 1
fi

if ! [[ "$MIN_DEPTH" =~ ^[0-9]+$ ]]; then
    echo "Error: Minimum depth must be a positive integer"
    exit 1
fi

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_VCF")
mkdir -p "$OUTPUT_DIR"

# Log file
LOG_FILE="${OUTPUT_VCF%.vcf.gz}.filter.log"
if [[ "$OUTPUT_VCF" == *.vcf ]]; then
    LOG_FILE="${OUTPUT_VCF%.vcf}.filter.log"
fi

echo "========================================" | tee "$LOG_FILE"
echo "Variant Filtering" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "Input VCF: $INPUT_VCF" | tee -a "$LOG_FILE"
echo "Output VCF: $OUTPUT_VCF" | tee -a "$LOG_FILE"
echo "Minimum QUAL: $MIN_QUAL" | tee -a "$LOG_FILE"
echo "Minimum depth: $MIN_DEPTH" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Activate conda environment
echo "Activating conda environment: nanopore_pipeline" | tee -a "$LOG_FILE"
eval "$(conda shell.bash hook)"
conda activate nanopore_pipeline

# Check if bcftools is available
if ! command -v bcftools &> /dev/null; then
    echo "Error: bcftools not found"
    exit 1
fi

# Count variants before filtering
echo "Counting variants before filtering..." | tee -a "$LOG_FILE"
BEFORE_COUNT=$(bcftools view -H "$INPUT_VCF" | wc -l)
echo "Variants before filtering: $BEFORE_COUNT" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Determine if output should be compressed
if [[ "$OUTPUT_VCF" == *.gz ]]; then
    OUTPUT_UNCOMPRESSED="${OUTPUT_VCF%.gz}"
    COMPRESS_OUTPUT=true
else
    OUTPUT_UNCOMPRESSED="$OUTPUT_VCF"
    COMPRESS_OUTPUT=false
fi

# Check if VCF has INFO/DP or FORMAT/DP
echo "Detecting depth field location..." | tee -a "$LOG_FILE"
HAS_INFO_DP=$(bcftools view -h "$INPUT_VCF" | grep -c "^##INFO=<ID=DP" || true)
HAS_FORMAT_DP=$(bcftools view -h "$INPUT_VCF" | grep -c "^##FORMAT=<ID=DP" || true)

echo "INFO/DP present: $([ $HAS_INFO_DP -gt 0 ] && echo 'Yes' || echo 'No')" | tee -a "$LOG_FILE"
echo "FORMAT/DP present: $([ $HAS_FORMAT_DP -gt 0 ] && echo 'Yes' || echo 'No')" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Apply filters based on available depth fields
echo "Applying filters..." | tee -a "$LOG_FILE"

if [ $HAS_INFO_DP -gt 0 ]; then
    # VCF has INFO/DP (e.g., from some callers)
    echo "Using INFO/DP for depth filtering" | tee -a "$LOG_FILE"
    FILTER_EXPR="QUAL>=${MIN_QUAL} && INFO/DP>=${MIN_DEPTH}"
elif [ $HAS_FORMAT_DP -gt 0 ]; then
    # VCF has FORMAT/DP (e.g., Clair3, many other callers)
    echo "Using FORMAT/DP for depth filtering" | tee -a "$LOG_FILE"
    FILTER_EXPR="QUAL>=${MIN_QUAL} && FORMAT/DP>=${MIN_DEPTH}"
else
    # No depth field found, filter on QUAL only
    echo "Warning: No DP field found, filtering on QUAL only" | tee -a "$LOG_FILE"
    FILTER_EXPR="QUAL>=${MIN_QUAL}"
fi

echo "Filter expression: $FILTER_EXPR" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Apply filter
bcftools view \
    -i "$FILTER_EXPR" \
    -O v \
    "$INPUT_VCF" > "$OUTPUT_UNCOMPRESSED" 2>> "$LOG_FILE"

if [ $? -ne 0 ]; then
    echo "Error: Filtering failed" | tee -a "$LOG_FILE"
    rm -f "$OUTPUT_UNCOMPRESSED"
    exit 1
fi

# Count variants after filtering
AFTER_COUNT=$(bcftools view -H "$OUTPUT_UNCOMPRESSED" | wc -l)
FILTERED_COUNT=$((BEFORE_COUNT - AFTER_COUNT))

if [ $BEFORE_COUNT -gt 0 ]; then
    PERCENT_KEPT=$(awk "BEGIN {printf \"%.2f\", ($AFTER_COUNT/$BEFORE_COUNT)*100}")
else
    PERCENT_KEPT="0.00"
fi

echo "" | tee -a "$LOG_FILE"
echo "Filtering completed!" | tee -a "$LOG_FILE"
echo "Variants before filtering: $BEFORE_COUNT" | tee -a "$LOG_FILE"
echo "Variants after filtering: $AFTER_COUNT" | tee -a "$LOG_FILE"
echo "Variants filtered out: $FILTERED_COUNT" | tee -a "$LOG_FILE"
echo "Percentage kept: ${PERCENT_KEPT}%" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Compress and index if requested
if [ "$COMPRESS_OUTPUT" = true ]; then
    echo "Compressing and indexing VCF..." | tee -a "$LOG_FILE"
    bgzip -f "$OUTPUT_UNCOMPRESSED"
    tabix -f -p vcf "$OUTPUT_VCF"
    echo "Output: $OUTPUT_VCF" | tee -a "$LOG_FILE"
    echo "Index: ${OUTPUT_VCF}.tbi" | tee -a "$LOG_FILE"
else
    echo "Output: $OUTPUT_VCF" | tee -a "$LOG_FILE"
fi

# Generate summary statistics
echo "" | tee -a "$LOG_FILE"
echo "Summary statistics for filtered variants:" | tee -a "$LOG_FILE"
echo "------------------------------------------" | tee -a "$LOG_FILE"

if [ "$COMPRESS_OUTPUT" = true ]; then
    VCF_TO_ANALYZE="$OUTPUT_VCF"
else
    VCF_TO_ANALYZE="$OUTPUT_VCF"
fi

bcftools stats "$VCF_TO_ANALYZE" | grep "^SN" | cut -f 3- | tee -a "$LOG_FILE"

echo "" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

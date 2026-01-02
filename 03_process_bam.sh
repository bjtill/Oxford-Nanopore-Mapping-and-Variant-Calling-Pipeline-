#!/bin/bash

################################################################################
# Script: 03_process_bam.sh
# Description: Sort, index, and generate statistics for BAM file
# Usage: ./03_process_bam.sh <input_bam> <threads>
################################################################################

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <input_bam> <threads>"
    echo ""
    echo "Arguments:"
    echo "  input_bam    Path to input BAM file"
    echo "  threads      Number of threads to use"
    echo ""
    echo "Output:"
    echo "  Creates <input_bam>.sorted.bam (sorted BAM)"
    echo "  Creates <input_bam>.sorted.bam.bai (BAM index)"
    echo "  Creates <input_bam>.sorted.stats.txt (alignment statistics)"
    echo ""
    echo "Example:"
    echo "  $0 aligned.bam 16"
    exit 1
}

# Check arguments
if [ $# -ne 2 ]; then
    usage
fi

INPUT_BAM="$1"
THREADS="$2"

# Validate inputs
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file not found: $INPUT_BAM"
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error: Threads must be a positive integer"
    exit 1
fi

# Output files
SORTED_BAM="${INPUT_BAM%.bam}.sorted.bam"
STATS_FILE="${SORTED_BAM%.bam}.stats.txt"
LOG_FILE="${SORTED_BAM%.bam}.process.log"

echo "========================================" | tee "$LOG_FILE"
echo "BAM Processing (Sort, Index, Stats)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "Input BAM: $INPUT_BAM" | tee -a "$LOG_FILE"
echo "Output sorted BAM: $SORTED_BAM" | tee -a "$LOG_FILE"
echo "Threads: $THREADS" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Activate conda environment
echo "Activating conda environment: nanopore_pipeline" | tee -a "$LOG_FILE"
eval "$(conda shell.bash hook)"
conda activate nanopore_pipeline

# Check if samtools is available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found"
    exit 1
fi

# Sort BAM file
echo "Sorting BAM file..." | tee -a "$LOG_FILE"
samtools sort \
    -@ "$THREADS" \
    -m 4G \
    -o "$SORTED_BAM" \
    "$INPUT_BAM" 2>&1 | tee -a "$LOG_FILE"

if [ $? -ne 0 ]; then
    echo "Error: BAM sorting failed" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Sorting completed successfully!" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Index BAM file
echo "Indexing BAM file..." | tee -a "$LOG_FILE"
samtools index \
    -@ "$THREADS" \
    "$SORTED_BAM" 2>&1 | tee -a "$LOG_FILE"

if [ $? -ne 0 ]; then
    echo "Error: BAM indexing failed" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Indexing completed successfully!" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Generate detailed statistics
echo "Generating alignment statistics..." | tee -a "$LOG_FILE"

{
    echo "=========================================="
    echo "Alignment Statistics"
    echo "=========================================="
    echo ""
    echo "Basic flagstat:"
    echo "------------------------------------------"
    samtools flagstat "$SORTED_BAM"
    echo ""
    echo "Detailed stats:"
    echo "------------------------------------------"
    samtools stats "$SORTED_BAM" | grep ^SN | cut -f 2-
    echo ""
    echo "Coverage statistics by chromosome:"
    echo "------------------------------------------"
    samtools idxstats "$SORTED_BAM"
    echo ""
    echo "Read length distribution:"
    echo "------------------------------------------"
    samtools view "$SORTED_BAM" | \
        awk '{print length($10)}' | \
        sort -n | \
        awk '{
            count[NR] = $1;
            sum += $1;
        }
        END {
            print "Mean read length: " sum/NR;
            print "Median read length: " count[int(NR/2)];
            print "Min read length: " count[1];
            print "Max read length: " count[NR];
        }'
} > "$STATS_FILE" 2>&1

echo "Statistics saved to: $STATS_FILE" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Display summary
echo "========================================" | tee -a "$LOG_FILE"
echo "Processing Summary" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Sorted BAM: $SORTED_BAM" | tee -a "$LOG_FILE"
echo "BAM index: ${SORTED_BAM}.bai" | tee -a "$LOG_FILE"
echo "Statistics: $STATS_FILE" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Show quick summary
echo "Quick summary:" | tee -a "$LOG_FILE"
samtools flagstat "$SORTED_BAM" | head -5 | tee -a "$LOG_FILE"

echo "" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

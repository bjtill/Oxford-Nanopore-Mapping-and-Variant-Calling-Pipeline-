#!/bin/bash

################################################################################
# Script: 04a_call_variants_clair3_triploid.sh
# Description: Integrated pipeline for triploid variant calling using Clair3 + FreeBayes
# Usage: ./04a_call_variants_clair3_triploid.sh <input_bam> <reference> <output_dir> <model> <threads> <ploidy> [min_qual] [min_depth]
################################################################################

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <input_bam> <reference_fasta> <output_dir> <model> <threads> <ploidy> [min_qual] [min_depth]"
    echo ""
    echo "Required Arguments:"
    echo "  input_bam        Path to sorted and indexed BAM file"
    echo "  reference_fasta  Path to reference genome FASTA file"
    echo "  output_dir       Directory for output"
    echo "  model            Clair3 model name:"
    echo "                   - r1041_e82_400bps_sup_v500 (R10.4.1, recommended)"
    echo "                   - r1041_e82_400bps_sup_v410 (R10.4.1, older)"
    echo "                   - r941_prom_sup_g5014 (R9.4.1)"
    echo "                   - r941_prom_hac_g360+g422 (R9.4.1)"
    echo "  threads          Number of threads to use"
    echo "  ploidy           Sample ploidy (2 for diploid, 3 for triploid)"
    echo ""
    echo "Optional Arguments:"
    echo "  min_qual         Minimum variant quality score (default: 20)"
    echo "  min_depth        Minimum read depth (default: 5)"
    echo ""
    echo "Examples:"
    echo "  # Triploid with default filters (QUAL>=20, DP>=5)"
    echo "  $0 sample.sorted.bam reference.fasta output r1041_e82_400bps_sup_v500 16 3"
    echo ""
    echo "  # Triploid with custom filters (QUAL>=30, DP>=10)"
    echo "  $0 sample.sorted.bam reference.fasta output r1041_e82_400bps_sup_v500 16 3 30 10"
    echo ""
    echo "  # Diploid with relaxed filters (QUAL>=15, DP>=3)"
    echo "  $0 sample.sorted.bam reference.fasta output r1041_e82_400bps_sup_v500 16 2 15 3"
    echo ""
    echo "Pipeline steps:"
    echo "  1. Clair3 variant calling (discovers variants)"
    echo "  2. Filter variants by QUAL and DP"
    echo "  3. For ploidy > 2: FreeBayes re-genotyping with correct ploidy"
    exit 1
}

# Check arguments
if [ $# -lt 6 ] || [ $# -gt 8 ]; then
    usage
fi

INPUT_BAM="$1"
REFERENCE="$2"
OUTPUT_DIR="$3"
MODEL="$4"
THREADS="$5"
PLOIDY="$6"

# Set filter defaults or use provided values
MIN_QUAL="${7:-20}"
MIN_DEPTH="${8:-5}"

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

if ! [[ "$PLOIDY" =~ ^[0-9]+$ ]] || [ "$PLOIDY" -lt 2 ]; then
    echo "Error: Ploidy must be an integer >= 2"
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

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Define subdirectories
CLAIR3_DIR="${OUTPUT_DIR}/01_clair3_discovery"
FILTERED_DIR="${OUTPUT_DIR}/02_filtered_variants"
FINAL_DIR="${OUTPUT_DIR}/03_final_genotypes"

# Log file
LOG_FILE="${OUTPUT_DIR}/triploid_pipeline.log"

echo "========================================" | tee "$LOG_FILE"
echo "Integrated Clair3 + FreeBayes Pipeline" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "Input BAM: $INPUT_BAM" | tee -a "$LOG_FILE"
echo "Reference: $REFERENCE" | tee -a "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR" | tee -a "$LOG_FILE"
echo "Model: $MODEL" | tee -a "$LOG_FILE"
echo "Threads: $THREADS" | tee -a "$LOG_FILE"
echo "Ploidy: $PLOIDY" | tee -a "$LOG_FILE"
echo "Min QUAL: $MIN_QUAL" | tee -a "$LOG_FILE"
echo "Min Depth: $MIN_DEPTH" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

###############################################################################
# STEP 1: Clair3 Variant Discovery
###############################################################################

echo "========================================" | tee -a "$LOG_FILE"
echo "STEP 1: Variant Discovery with Clair3" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Activate Clair3 environment
echo "Activating conda environment: clair3" | tee -a "$LOG_FILE"
eval "$(conda shell.bash hook)"
conda activate clair3

# Check if Clair3 is available
if ! command -v run_clair3.sh &> /dev/null; then
    echo "Error: Clair3 not found. Please install it in the clair3 conda environment."
    exit 1
fi

# Get model path
CLAIR3_BIN=$(dirname $(which run_clair3.sh))
MODEL_PATH="${CLAIR3_BIN}/models/${MODEL}"

if [ ! -d "$MODEL_PATH" ]; then
    echo "Error: Model not found at: $MODEL_PATH" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    echo "Available models in ${CLAIR3_BIN}/models:" | tee -a "$LOG_FILE"
    ls -1 "${CLAIR3_BIN}/models/" 2>/dev/null | grep -v "^ont" | tee -a "$LOG_FILE"
    exit 1
fi

echo "Using model: $MODEL_PATH" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Run Clair3
if [ ! -f "${CLAIR3_DIR}/merge_output.vcf.gz" ]; then
    echo "Running Clair3 variant calling..." | tee -a "$LOG_FILE"
    
    run_clair3.sh \
        --bam_fn="$INPUT_BAM" \
        --ref_fn="$REFERENCE" \
        --output="$CLAIR3_DIR" \
        --threads="$THREADS" \
        --platform=ont \
        --model_path="$MODEL_PATH" \
        --include_all_ctgs \
        --no_phasing_for_fa 2>&1 | tee -a "$LOG_FILE"
    
    if [ $? -ne 0 ]; then
        echo "Error: Clair3 failed" | tee -a "$LOG_FILE"
        exit 1
    fi
    
    echo "Clair3 completed successfully!" | tee -a "$LOG_FILE"
else
    echo "Clair3 output already exists, skipping..." | tee -a "$LOG_FILE"
fi

# Quick variant count
CLAIR3_COUNT=$(zgrep -v "^#" "${CLAIR3_DIR}/merge_output.vcf.gz" | wc -l)
echo "Variants discovered by Clair3: $CLAIR3_COUNT" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

###############################################################################
# STEP 2: Filter Variants
###############################################################################

echo "========================================" | tee -a "$LOG_FILE"
echo "STEP 2: Filter High-Quality Variants" | tee -a "$LOG_FILE"
echo "       (QUAL>=${MIN_QUAL}, DP>=${MIN_DEPTH})" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Switch to nanopore_pipeline environment for bcftools
conda activate nanopore_pipeline

mkdir -p "$FILTERED_DIR"
FILTERED_VCF="${FILTERED_DIR}/clair3_filtered.vcf.gz"

if [ ! -f "$FILTERED_VCF" ]; then
    echo "Filtering variants (QUAL>=${MIN_QUAL}, DP>=${MIN_DEPTH})..." | tee -a "$LOG_FILE"
    
    bcftools view \
        -i "QUAL>=${MIN_QUAL} && FORMAT/DP>=${MIN_DEPTH}" \
        -O z \
        -o "$FILTERED_VCF" \
        "${CLAIR3_DIR}/merge_output.vcf.gz" 2>&1 | tee -a "$LOG_FILE"
    
    tabix -f -p vcf "$FILTERED_VCF"
    
    echo "Filtering completed!" | tee -a "$LOG_FILE"
else
    echo "Filtered variants already exist, skipping..." | tee -a "$LOG_FILE"
fi

FILTERED_COUNT=$(zgrep -v "^#" "$FILTERED_VCF" | wc -l)
PERCENT_KEPT=$(awk "BEGIN {printf \"%.2f\", ($FILTERED_COUNT/$CLAIR3_COUNT)*100}")
echo "Variants after filtering: $FILTERED_COUNT (${PERCENT_KEPT}% of discovered)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

###############################################################################
# STEP 3: Re-genotype with FreeBayes (if ploidy > 2)
###############################################################################

if [ "$PLOIDY" -eq 2 ]; then
    echo "========================================" | tee -a "$LOG_FILE"
    echo "Ploidy = 2 (diploid), skipping FreeBayes" | tee -a "$LOG_FILE"
    echo "========================================" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    echo "Final VCF: $FILTERED_VCF" | tee -a "$LOG_FILE"
    
    # Create symlink to final output
    mkdir -p "$FINAL_DIR"
    ln -sf "../02_filtered_variants/clair3_filtered.vcf.gz" "${FINAL_DIR}/final_diploid.vcf.gz"
    ln -sf "../02_filtered_variants/clair3_filtered.vcf.gz.tbi" "${FINAL_DIR}/final_diploid.vcf.gz.tbi"
    
else
    echo "========================================" | tee -a "$LOG_FILE"
    echo "STEP 3: Re-genotype with FreeBayes" | tee -a "$LOG_FILE"
    echo "       (Ploidy = $PLOIDY)" | tee -a "$LOG_FILE"
    echo "========================================" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    
    # Check if FreeBayes is available
    if ! command -v freebayes &> /dev/null; then
        echo "Error: FreeBayes not found. Installing..." | tee -a "$LOG_FILE"
        mamba install -c bioconda freebayes -y 2>&1 | tee -a "$LOG_FILE"
    fi
    
    mkdir -p "$FINAL_DIR"
    
    # Extract variant sites for FreeBayes
    SITES_FILE="${FILTERED_DIR}/variant_sites.bed"
    if [ ! -f "$SITES_FILE" ]; then
        echo "Extracting variant sites..." | tee -a "$LOG_FILE"
        bcftools query -f '%CHROM\t%POS0\t%POS\n' "$FILTERED_VCF" > "$SITES_FILE"
        echo "Extracted $(wc -l < $SITES_FILE) variant sites" | tee -a "$LOG_FILE"
    fi
    
    # Run FreeBayes
    FREEBAYES_VCF="${FINAL_DIR}/freebayes_ploidy${PLOIDY}.vcf"
    
    if [ ! -f "${FREEBAYES_VCF}.gz" ]; then
        echo "Running FreeBayes with ploidy=$PLOIDY..." | tee -a "$LOG_FILE"
        echo "This may take several hours..." | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
        
        # Calculate min alternate fraction based on ploidy
        # For triploid: 1/3 ≈ 0.33, use 0.25 with buffer
        MIN_ALT_FRAC=$(awk "BEGIN {printf \"%.3f\", (1/$PLOIDY) - 0.08}")
        
        echo "FreeBayes parameters:" | tee -a "$LOG_FILE"
        echo "  Ploidy: $PLOIDY" | tee -a "$LOG_FILE"
        echo "  Min alternate fraction: $MIN_ALT_FRAC" | tee -a "$LOG_FILE"
        echo "  Min coverage: 10" | tee -a "$LOG_FILE"
        echo "" | tee -a "$LOG_FILE"
        
        freebayes \
            --fasta-reference "$REFERENCE" \
            --ploidy "$PLOIDY" \
            --min-alternate-fraction "$MIN_ALT_FRAC" \
            --min-coverage 10 \
            --min-mapping-quality 20 \
            --min-base-quality 10 \
            --targets "$SITES_FILE" \
            --genotype-qualities \
            "$INPUT_BAM" > "$FREEBAYES_VCF" 2>> "$LOG_FILE"
        
        if [ $? -ne 0 ]; then
            echo "Error: FreeBayes failed" | tee -a "$LOG_FILE"
            exit 1
        fi
        
        # Compress and index
        bgzip -f "$FREEBAYES_VCF"
        tabix -f -p vcf "${FREEBAYES_VCF}.gz"
        
        echo "FreeBayes completed successfully!" | tee -a "$LOG_FILE"
    else
        echo "FreeBayes output already exists, skipping..." | tee -a "$LOG_FILE"
    fi
    
    FREEBAYES_COUNT=$(zgrep -v "^#" "${FREEBAYES_VCF}.gz" | wc -l)
    FB_PERCENT=$(awk "BEGIN {printf \"%.2f\", ($FREEBAYES_COUNT/$FILTERED_COUNT)*100}")
    echo "Variants called by FreeBayes: $FREEBAYES_COUNT (${FB_PERCENT}% of filtered)" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
    
    # Create final symlink
    ln -sf "freebayes_ploidy${PLOIDY}.vcf.gz" "${FINAL_DIR}/final_ploidy${PLOIDY}.vcf.gz"
    ln -sf "freebayes_ploidy${PLOIDY}.vcf.gz.tbi" "${FINAL_DIR}/final_ploidy${PLOIDY}.vcf.gz.tbi"
fi

###############################################################################
# STEP 4: Generate Summary Statistics
###############################################################################

echo "========================================" | tee -a "$LOG_FILE"
echo "STEP 4: Summary Statistics" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

SUMMARY_FILE="${OUTPUT_DIR}/pipeline_summary.txt"

{
    echo "=========================================="
    echo "VARIANT CALLING PIPELINE SUMMARY"
    echo "=========================================="
    echo "Date: $(date)"
    echo ""
    echo "Sample: $INPUT_BAM"
    echo "Reference: $REFERENCE"
    echo "Ploidy: $PLOIDY"
    echo ""
    echo "=========================================="
    echo "FILTER PARAMETERS"
    echo "=========================================="
    echo ""
    echo "Minimum QUAL: $MIN_QUAL"
    echo "Minimum Depth: $MIN_DEPTH"
    echo ""
    echo "=========================================="
    echo "VARIANT COUNTS"
    echo "=========================================="
    echo ""
    echo "1. Clair3 Discovery:    $CLAIR3_COUNT variants"
    echo "2. After Filtering:     $FILTERED_COUNT variants (${PERCENT_KEPT}%)"
    
    if [ "$PLOIDY" -gt 2 ]; then
        echo "3. FreeBayes (ploidy=$PLOIDY): $FREEBAYES_COUNT variants (${FB_PERCENT}% of filtered)"
    fi
    
    echo ""
    echo "=========================================="
    echo "OUTPUT FILES"
    echo "=========================================="
    echo ""
    echo "Clair3 VCF:     ${CLAIR3_DIR}/merge_output.vcf.gz"
    echo "Filtered VCF:   $FILTERED_VCF"
    
    if [ "$PLOIDY" -eq 2 ]; then
        echo "Final VCF:      ${FINAL_DIR}/final_diploid.vcf.gz"
    else
        echo "FreeBayes VCF:  ${FREEBAYES_VCF}.gz"
        echo "Final VCF:      ${FINAL_DIR}/final_ploidy${PLOIDY}.vcf.gz"
    fi
    
    echo ""
    echo "=========================================="
    echo "GENOTYPE INTERPRETATION"
    echo "=========================================="
    echo ""
    
    if [ "$PLOIDY" -eq 2 ]; then
        echo "Diploid genotypes:"
        echo "  0/0 - Homozygous reference"
        echo "  0/1 - Heterozygous"
        echo "  1/1 - Homozygous alternate"
    elif [ "$PLOIDY" -eq 3 ]; then
        echo "Triploid genotypes:"
        echo "  0/0/0 - Homozygous reference"
        echo "  0/0/1 - One alternate allele (AAB)"
        echo "  0/1/1 - Two alternate alleles (ABB)"
        echo "  1/1/1 - Homozygous alternate"
    else
        echo "Ploidy $PLOIDY genotypes will have $PLOIDY alleles"
    fi
    
    echo ""
    
} > "$SUMMARY_FILE"

cat "$SUMMARY_FILE" | tee -a "$LOG_FILE"

echo "" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Pipeline Completed Successfully!" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Summary: $SUMMARY_FILE" | tee -a "$LOG_FILE"

if [ "$PLOIDY" -eq 2 ]; then
    echo "Final VCF: ${FINAL_DIR}/final_diploid.vcf.gz" | tee -a "$LOG_FILE"
else
    echo "Final VCF: ${FINAL_DIR}/final_ploidy${PLOIDY}.vcf.gz" | tee -a "$LOG_FILE"
fi

echo "========================================" | tee -a "$LOG_FILE"

#!/bin/bash

################################################################################
# Script: 06_variant_stats.sh
# Description: Generate comprehensive statistics for variant calls
# Usage: ./06_variant_stats.sh <input_vcf> <reference> <output_prefix>
################################################################################

set -euo pipefail

# Function to display usage
usage() {
    echo "Usage: $0 <input_vcf> <reference_fasta> <output_prefix>"
    echo ""
    echo "Arguments:"
    echo "  input_vcf        Path to VCF file (can be gzipped)"
    echo "  reference_fasta  Path to reference genome FASTA file"
    echo "  output_prefix    Prefix for output statistics files"
    echo ""
    echo "Example:"
    echo "  $0 variants.vcf.gz reference.fasta variant_stats"
    echo ""
    echo "Output files:"
    echo "  <prefix>.summary.txt      - Overall summary statistics"
    echo "  <prefix>.by_type.txt      - Variant counts by type"
    echo "  <prefix>.by_chromosome.txt - Variant distribution by chromosome"
    echo "  <prefix>.quality_dist.txt - Quality score distribution"
    echo "  <prefix>.depth_dist.txt   - Depth distribution"
    exit 1
}

# Check arguments
if [ $# -ne 3 ]; then
    usage
fi

INPUT_VCF="$1"
REFERENCE="$2"
OUTPUT_PREFIX="$3"

# Validate inputs
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input VCF file not found: $INPUT_VCF"
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference FASTA file not found: $REFERENCE"
    exit 1
fi

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
mkdir -p "$OUTPUT_DIR"

# Output files
SUMMARY_FILE="${OUTPUT_PREFIX}.summary.txt"
BY_TYPE_FILE="${OUTPUT_PREFIX}.by_type.txt"
BY_CHR_FILE="${OUTPUT_PREFIX}.by_chromosome.txt"
QUAL_DIST_FILE="${OUTPUT_PREFIX}.quality_dist.txt"
DEPTH_DIST_FILE="${OUTPUT_PREFIX}.depth_dist.txt"
LOG_FILE="${OUTPUT_PREFIX}.stats.log"

echo "========================================" | tee "$LOG_FILE"
echo "Variant Statistics" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "Input VCF: $INPUT_VCF" | tee -a "$LOG_FILE"
echo "Reference: $REFERENCE" | tee -a "$LOG_FILE"
echo "Output prefix: $OUTPUT_PREFIX" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Activate conda environment
echo "Activating conda environment: nanopore_pipeline" | tee -a "$LOG_FILE"
eval "$(conda shell.bash hook)"
conda activate nanopore_pipeline

# Check if required tools are available
if ! command -v bcftools &> /dev/null; then
    echo "Error: bcftools not found"
    exit 1
fi

if ! command -v rtg &> /dev/null; then
    echo "Warning: rtg-tools not found. Some statistics will be skipped." | tee -a "$LOG_FILE"
    RTG_AVAILABLE=false
else
    RTG_AVAILABLE=true
fi

# Count different variant types (for use throughout script)
echo "Counting variant types..." | tee -a "$LOG_FILE"
TOTAL_VARIANTS=$(bcftools view -H "$INPUT_VCF" 2>/dev/null | wc -l)
SNP_COUNT=$(bcftools view -v snps "$INPUT_VCF" 2>/dev/null | bcftools query -f '%CHROM\n' | wc -l)
INSERTION_COUNT=$(bcftools view -v indels "$INPUT_VCF" 2>/dev/null | bcftools query -i 'strlen(ALT)>strlen(REF)' -f '%CHROM\n' | wc -l)
DELETION_COUNT=$(bcftools view -v indels "$INPUT_VCF" 2>/dev/null | bcftools query -i 'strlen(ALT)<strlen(REF)' -f '%CHROM\n' | wc -l)
MNP_COUNT=$(bcftools view -v mnps "$INPUT_VCF" 2>/dev/null | bcftools query -f '%CHROM\n' | wc -l)
OTHER_COUNT=$(bcftools view -v other "$INPUT_VCF" 2>/dev/null | bcftools query -f '%CHROM\n' | wc -l)

echo "Total variants: $TOTAL_VARIANTS" | tee -a "$LOG_FILE"
echo "  SNPs: $SNP_COUNT" | tee -a "$LOG_FILE"
echo "  Insertions: $INSERTION_COUNT" | tee -a "$LOG_FILE"
echo "  Deletions: $DELETION_COUNT" | tee -a "$LOG_FILE"
echo "  MNPs: $MNP_COUNT" | tee -a "$LOG_FILE"
echo "  Other: $OTHER_COUNT" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Generate overall summary statistics
echo "Generating overall summary statistics..." | tee -a "$LOG_FILE"

{
    echo "=========================================="
    echo "OVERALL VARIANT SUMMARY"
    echo "=========================================="
    echo "Analysis date: $(date)"
    echo "Input VCF: $INPUT_VCF"
    echo ""
    
    # Basic counts
    echo "BASIC COUNTS:"
    echo "-------------"
    bcftools stats "$INPUT_VCF" 2>/dev/null | grep "^SN" | cut -f 3-
    echo ""
    
    # Variant types with descriptions
    echo "VARIANT TYPES:"
    echo "-------------"
    echo "SNPs (Single Nucleotide Polymorphisms): $SNP_COUNT"
    echo "Insertions: $INSERTION_COUNT"
    echo "Deletions: $DELETION_COUNT"
    echo "MNPs (Multi-Nucleotide Polymorphisms): $MNP_COUNT"
    echo "Other (complex/symbolic variants): $OTHER_COUNT"
    echo ""
    echo "Total indels (insertions + deletions): $((INSERTION_COUNT + DELETION_COUNT))"
    echo ""
    
    # Transition/Transversion ratio
    echo "Ti/Tv RATIO:"
    echo "-----------"
    bcftools stats "$INPUT_VCF" 2>/dev/null | grep "^TSTV" | cut -f 2-
    echo ""
    
    # Quality statistics
    echo "QUALITY STATISTICS:"
    echo "------------------"
    bcftools query -f '%QUAL\n' "$INPUT_VCF" | awk '
    BEGIN {sum=0; count=0; min=999999; max=0}
    {
        sum+=$1; count++;
        if($1<min) min=$1;
        if($1>max) max=$1;
        qual[count]=$1;
    }
    END {
        asort(qual);
        median = (count%2==0) ? (qual[count/2]+qual[count/2+1])/2 : qual[int(count/2)+1];
        printf "Mean QUAL: %.2f\n", sum/count;
        printf "Median QUAL: %.2f\n", median;
        printf "Min QUAL: %.2f\n", min;
        printf "Max QUAL: %.2f\n", max;
    }'
    echo ""
    
    # Depth statistics (if available)
    if bcftools query -f '%INFO/DP\n' "$INPUT_VCF" 2>/dev/null | head -1 | grep -q "[0-9]"; then
        echo "DEPTH STATISTICS:"
        echo "----------------"
        bcftools query -f '%INFO/DP\n' "$INPUT_VCF" | awk '
        BEGIN {sum=0; count=0; min=999999; max=0}
        {
            sum+=$1; count++;
            if($1<min) min=$1;
            if($1>max) max=$1;
            dp[count]=$1;
        }
        END {
            asort(dp);
            median = (count%2==0) ? (dp[count/2]+dp[count/2+1])/2 : dp[int(count/2)+1];
            printf "Mean depth: %.2f\n", sum/count;
            printf "Median depth: %.2f\n", median;
            printf "Min depth: %d\n", min;
            printf "Max depth: %d\n", max;
        }'
        echo ""
    fi
} > "$SUMMARY_FILE" 2>&1

echo "Summary statistics saved to: $SUMMARY_FILE" | tee -a "$LOG_FILE"

# Generate variant counts by type (with descriptions)
echo "Generating variant type distribution..." | tee -a "$LOG_FILE"

TEMP_BY_TYPE=$(mktemp)
{
    echo "VARIANT TYPE DISTRIBUTION"
    echo "========================="
    echo ""
    echo "Variant Type Definitions:"
    echo "  SNP       = Single Nucleotide Polymorphism (1 base change)"
    echo "  Insertion = Added bases relative to reference"
    echo "  Deletion  = Removed bases relative to reference"
    echo "  MNP       = Multi-Nucleotide Polymorphism (multiple consecutive base changes)"
    echo "  Other     = Complex or symbolic variants (breakends, large SVs, etc.)"
    echo ""
    echo -e "Type\tCount\tPercentage"
    echo -e "SNP\t${SNP_COUNT}\t$(awk "BEGIN {printf \"%.2f\", ($SNP_COUNT/$TOTAL_VARIANTS)*100}")%"
    echo -e "Insertion\t${INSERTION_COUNT}\t$(awk "BEGIN {printf \"%.2f\", ($INSERTION_COUNT/$TOTAL_VARIANTS)*100}")%"
    echo -e "Deletion\t${DELETION_COUNT}\t$(awk "BEGIN {printf \"%.2f\", ($DELETION_COUNT/$TOTAL_VARIANTS)*100}")%"
    echo -e "MNP\t${MNP_COUNT}\t$(awk "BEGIN {printf \"%.2f\", ($MNP_COUNT/$TOTAL_VARIANTS)*100}")%"
    echo -e "Other\t${OTHER_COUNT}\t$(awk "BEGIN {printf \"%.2f\", ($OTHER_COUNT/$TOTAL_VARIANTS)*100}")%"
    echo -e "---\t---\t---"
    echo -e "TOTAL\t${TOTAL_VARIANTS}\t100.00%"
} > "$TEMP_BY_TYPE"
mv "$TEMP_BY_TYPE" "$BY_TYPE_FILE"

echo "Variant type distribution saved to: $BY_TYPE_FILE" | tee -a "$LOG_FILE"

# Generate variant distribution by chromosome (with all variant types)
echo "Generating chromosome distribution..." | tee -a "$LOG_FILE"

TEMP_BY_CHR=$(mktemp)
{
    echo "VARIANT DISTRIBUTION BY CHROMOSOME"
    echo "=================================="
    echo ""
    echo -e "Chromosome\tTotal\tSNPs\tInsertions\tDeletions\tMNPs\tOther"
    
    # Get all chromosomes
    bcftools query -f '%CHROM\n' "$INPUT_VCF" | sort -u | while read chr; do
        # Count total variants for this chromosome
        total=$(bcftools view -r "$chr" "$INPUT_VCF" 2>/dev/null | bcftools query -f '%CHROM\n' | wc -l)
        
        # Count SNPs
        snps=$(bcftools view -v snps -r "$chr" "$INPUT_VCF" 2>/dev/null | bcftools query -f '%CHROM\n' | wc -l)
        
        # Count insertions
        insertions=$(bcftools view -v indels -r "$chr" "$INPUT_VCF" 2>/dev/null | bcftools query -i 'strlen(ALT)>strlen(REF)' -f '%CHROM\n' | wc -l)
        
        # Count deletions
        deletions=$(bcftools view -v indels -r "$chr" "$INPUT_VCF" 2>/dev/null | bcftools query -i 'strlen(ALT)<strlen(REF)' -f '%CHROM\n' | wc -l)
        
        # Count MNPs
        mnps=$(bcftools view -v mnps -r "$chr" "$INPUT_VCF" 2>/dev/null | bcftools query -f '%CHROM\n' | wc -l)
        
        # Count Other
        other=$(bcftools view -v other -r "$chr" "$INPUT_VCF" 2>/dev/null | bcftools query -f '%CHROM\n' | wc -l)
        
        echo -e "${chr}\t${total}\t${snps}\t${insertions}\t${deletions}\t${mnps}\t${other}"
    done
    
    echo ""
    echo "Summary totals across all chromosomes:"
    echo -e "TOTAL\t${TOTAL_VARIANTS}\t${SNP_COUNT}\t${INSERTION_COUNT}\t${DELETION_COUNT}\t${MNP_COUNT}\t${OTHER_COUNT}"
} > "$TEMP_BY_CHR"
mv "$TEMP_BY_CHR" "$BY_CHR_FILE"

echo "Chromosome distribution saved to: $BY_CHR_FILE" | tee -a "$LOG_FILE"

# Generate quality score distribution
echo "Generating quality score distribution..." | tee -a "$LOG_FILE"

TEMP_QUAL=$(mktemp)
{
    echo "QUALITY SCORE DISTRIBUTION"
    echo "=========================="
    echo ""
    echo -e "QUAL_Range\tCount"
    bcftools query -f '%QUAL\n' "$INPUT_VCF" | awk '
    {
        bin = int($1/10)*10;
        count[bin]++;
    }
    END {
        for (b in count) {
            printf "%d-%d\t%d\n", b, b+9, count[b];
        }
    }' | sort -n
} > "$TEMP_QUAL"
mv "$TEMP_QUAL" "$QUAL_DIST_FILE"

echo "Quality distribution saved to: $QUAL_DIST_FILE" | tee -a "$LOG_FILE"

# Generate depth distribution
echo "Generating depth distribution..." | tee -a "$LOG_FILE"

if bcftools query -f '%INFO/DP\n' "$INPUT_VCF" 2>/dev/null | head -1 | grep -q "[0-9]"; then
    TEMP_DEPTH=$(mktemp)
    {
        echo "DEPTH DISTRIBUTION"
        echo "=================="
        echo ""
        echo -e "Depth_Range\tCount"
        bcftools query -f '%INFO/DP\n' "$INPUT_VCF" | awk '
        {
            if ($1 < 10) bin = int($1);
            else if ($1 < 50) bin = int($1/5)*5;
            else if ($1 < 100) bin = int($1/10)*10;
            else bin = int($1/50)*50;
            count[bin]++;
        }
        END {
            for (b in count) {
                if (b < 10) printf "%d\t%d\n", b, count[b];
                else if (b < 50) printf "%d-%d\t%d\n", b, b+4, count[b];
                else if (b < 100) printf "%d-%d\t%d\n", b, b+9, count[b];
                else printf "%d+\t%d\n", b, count[b];
            }
        }' | sort -n
    } > "$TEMP_DEPTH"
    mv "$TEMP_DEPTH" "$DEPTH_DIST_FILE"
    echo "Depth distribution saved to: $DEPTH_DIST_FILE" | tee -a "$LOG_FILE"
else
    echo "No depth information available in VCF" | tee -a "$LOG_FILE"
fi

# Run RTG vcfstats if available
if [ "$RTG_AVAILABLE" = true ]; then
    echo "" | tee -a "$LOG_FILE"
    echo "Running RTG vcfstats for additional metrics..." | tee -a "$LOG_FILE"
    rtg vcfstats "$INPUT_VCF" > "${OUTPUT_PREFIX}.rtg_stats.txt" 2>&1
    echo "RTG statistics saved to: ${OUTPUT_PREFIX}.rtg_stats.txt" | tee -a "$LOG_FILE"
fi

echo "" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"
echo "All statistics generated successfully!" | tee -a "$LOG_FILE"
echo "Output files:" | tee -a "$LOG_FILE"
echo "  - $SUMMARY_FILE" | tee -a "$LOG_FILE"
echo "  - $BY_TYPE_FILE" | tee -a "$LOG_FILE"
echo "  - $BY_CHR_FILE" | tee -a "$LOG_FILE"
echo "  - $QUAL_DIST_FILE" | tee -a "$LOG_FILE"
if [ -f "$DEPTH_DIST_FILE" ]; then
    echo "  - $DEPTH_DIST_FILE" | tee -a "$LOG_FILE"
fi
if [ -f "${OUTPUT_PREFIX}.rtg_stats.txt" ]; then
    echo "  - ${OUTPUT_PREFIX}.rtg_stats.txt" | tee -a "$LOG_FILE"
fi
echo "" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

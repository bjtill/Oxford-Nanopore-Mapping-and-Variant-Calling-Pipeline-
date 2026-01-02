# YAD GUI Integration Guide - Updated for v2

This guide provides complete information for integrating the Oxford Nanopore variant calling pipeline v2 scripts with a YAD (Yet Another Dialog) GUI interface.

NOTE: I hope to eventually add a YAD script for a GUI to run this pipeline.

## Overview of Pipeline v2

### Key Changes from v1:
1. **04a script now has TWO versions:**
   - `04a_call_variants_clair3_updated.sh` - Diploid only
   - `04a_call_variants_clair3_triploid_v2.sh` - Any ploidy (replaces 05 filtering)
2. **Filtering integrated** into triploid script
3. **05 script optional** for diploid workflow
4. **06 script enhanced** with complete variant type breakdown

### Recommended GUI Structure:

```
Main Menu:
├── [Sample Type Selection]
│   ├── Diploid (2n) → Use standard workflow
│   └── Polyploid (3n, 4n, etc.) → Use integrated workflow
├── [Step 1] Quality Control
├── [Step 2] Alignment
├── [Step 3] BAM Processing
├── [Step 4] Variant Calling
│   ├── [4a] Clair3 (diploid or polyploid)
│   ├── [4c] Longshot (SNVs only)
│   └── [Comparison] Compare multiple callers
├── [Step 5] Filtering (diploid only)
└── [Step 6] Statistics
```

---

## Script Parameters Reference

### 01_qc.sh - Quality Control
**Command:** `./01_qc.sh <input_fastq> <output_dir> <threads>`

**YAD Form Fields:**
```bash
--field="Input FASTQ:FL" "" \
--field="Output Directory:DIR" "" \
--field="Threads:NUM" "16!1..64!1"
```

**Example YAD Dialog:**
```bash
RESULT=$(yad --title="Step 1: Quality Control" \
    --form --width=600 \
    --field="Input FASTQ file:FL" "" \
    --field="Output directory:DIR" "" \
    --field="Number of threads:NUM" "16!1..64!1" \
    --button="Run:0" \
    --button="Cancel:1")
```

---

### 02_align.sh - Read Alignment
**Command:** `./02_align.sh <input_fastq> <reference> <output_bam> <threads>`

**YAD Form Fields:**
```bash
--field="Input FASTQ:FL" "" \
--field="Reference FASTA:FL" "" \
--field="Output BAM:SFL" "" \
--field="Threads:NUM" "16!1..64!1"
```

---

### 03_process_bam.sh - BAM Processing
**Command:** `./03_process_bam.sh <input_bam> <threads>`

**YAD Form Fields:**
```bash
--field="Input BAM:FL" "" \
--field="Threads:NUM" "16!1..64!1"
```

**Auto-Output Note:** Automatically creates `<input>.sorted.bam`

---

### 04a_call_variants_clair3_updated.sh - Diploid Variant Calling
**Command:** `./04a_call_variants_clair3_updated.sh <bam> <reference> <output_dir> <model> <threads>`

**YAD Form Fields:**
```bash
--field="Input BAM (sorted):FL" "" \
--field="Reference FASTA:FL" "" \
--field="Output Directory:DIR" "" \
--field="Clair3 Model:CB" "r1041_e82_400bps_sup_v500!r1041_e82_400bps_sup_v410!r941_prom_sup_g5014!r941_prom_hac_g360+g422" \
--field="Threads:NUM" "16!1..64!1"
```

**Model Selection Help Text:**
- r1041_e82_400bps_sup_v500: R10.4.1, super accuracy (recommended for R10)
- r1041_e82_400bps_sup_v410: R10.4.1, older version
- r941_prom_sup_g5014: R9.4.1, super accuracy
- r941_prom_hac_g360+g422: R9.4.1, high accuracy

**Use For:** Standard diploid samples (2n)

---

### 04a_call_variants_clair3_triploid_v2.sh - Polyploid Variant Calling (INTEGRATED)
**Command:** `./04a_call_variants_clair3_triploid_v2.sh <bam> <reference> <output_dir> <model> <threads> <ploidy> [min_qual] [min_depth]`

**YAD Form Fields:**
```bash
--field="Input BAM (sorted):FL" "" \
--field="Reference FASTA:FL" "" \
--field="Output Directory:DIR" "" \
--field="Clair3 Model:CB" "r1041_e82_400bps_sup_v500!r1041_e82_400bps_sup_v410!r941_prom_sup_g5014!r941_prom_hac_g360+g422" \
--field="Threads:NUM" "16!1..64!1" \
--field="Ploidy:NUM" "3!2..8!1" \
--field="Min QUAL (optional):NUM" "20!10..50!1" \
--field="Min Depth (optional):NUM" "5!3..20!1"
```

**Ploidy Options:**
- 2: Diploid (skips FreeBayes, uses Clair3 only)
- 3: Triploid (banana, seedless watermelon)
- 4: Tetraploid (potato, cotton)
- 6: Hexaploid (wheat)
- 8: Octoploid (strawberry)

**Filter Presets for GUI:**
```bash
--field="Filter Preset:CB" "Standard (QUAL≥20, DP≥5)!Relaxed (QUAL≥15, DP≥3)!Strict (QUAL≥30, DP≥10)!Publication (QUAL≥40, DP≥15)!Custom"
```

Then enable/disable manual filter fields based on selection.

**Use For:** Any non-diploid sample (3n, 4n, 6n, etc.)

**Important:** This script includes filtering - no need for separate filtering step!

---

### 04c_call_variants_longshot_updated.sh - Fast SNV Calling
**Command:** `./04c_call_variants_longshot_updated.sh <bam> <reference> <output_vcf> <threads>`

**YAD Form Fields:**
```bash
--field="Input BAM (sorted):FL" "" \
--field="Reference FASTA:FL" "" \
--field="Output VCF:SFL" "" \
--field="Threads:NUM" "8!1..32!1"
```

**Note:** Longshot is SNV-only, no indels

---

### 05_filter_variants_updated.sh - Variant Filtering (Diploid Workflow Only)
**Command:** `./05_filter_variants_updated.sh <input_vcf> <output_vcf> <min_qual> <min_depth>`

**YAD Form Fields:**
```bash
--field="Input VCF:FL" "" \
--field="Output VCF:SFL" "" \
--field="Min QUAL:NUM" "20!10..50!1" \
--field="Min Depth:NUM" "5!3..20!1"
```

**Filter Preset Dropdown:**
```bash
--field="Preset:CB" "Standard!Relaxed!Strict!Publication!Custom"
```

**Preset Values:**
- Standard: QUAL≥20, DP≥5
- Relaxed: QUAL≥15, DP≥3  
- Strict: QUAL≥30, DP≥10
- Publication: QUAL≥40, DP≥15
- Custom: User enters values

**Note:** Skip this step for polyploid samples (filtering is integrated into 04a_triploid)

---

### 06_variant_stats_v2.sh - Comprehensive Statistics
**Command:** `./06_variant_stats_v2.sh <input_vcf> <reference> <output_prefix>`

**YAD Form Fields:**
```bash
--field="Input VCF:FL" "" \
--field="Reference FASTA:FL" "" \
--field="Output Prefix:TXT" ""
```

**Output Files Created:**
- `<prefix>.summary.txt`
- `<prefix>.by_type.txt`
- `<prefix>.by_chromosome.txt`
- `<prefix>.quality_dist.txt`
- `<prefix>.depth_dist.txt`

---

### compare_vcfs.sh - Two-Way Comparison
**Command:** `./compare_vcfs.sh <vcf1> <vcf2> <output_prefix>`

**YAD Form Fields:**
```bash
--field="VCF 1 (e.g., Clair3):FL" "" \
--field="VCF 2 (e.g., Longshot):FL" "" \
--field="Output Prefix:TXT" ""
```

---

### compare_three_vcfs.sh - Three-Way Comparison
**Command:** `./compare_three_vcfs.sh <vcf1> <vcf2> <vcf3> <output_prefix>`

**YAD Form Fields:**
```bash
--field="VCF 1 (e.g., Clair3):FL" "" \
--field="VCF 2 (e.g., Longshot):FL" "" \
--field="VCF 3 (e.g., Pepper):FL" "" \
--field="Output Prefix:TXT" ""
```

---

## Recommended GUI Workflows

### Workflow 1: Diploid Sample (Standard)

```
┌─────────────────────────────────────┐
│   Diploid Workflow (2n)             │
├─────────────────────────────────────┤
│ 1. Quality Control (01_qc.sh)      │
│ 2. Alignment (02_align.sh)         │
│ 3. BAM Processing (03_process_bam) │
│ 4. Variant Calling (04a_diploid)   │
│ 5. Filter Variants (05_filter)     │
│ 6. Statistics (06_stats)           │
└─────────────────────────────────────┘
```

**YAD Implementation:**
```bash
# Main menu
yad --list --title="ONT Pipeline" \
    --column="Step" --column="Description" \
    "1" "Quality Control" \
    "2" "Alignment" \
    "3" "Process BAM" \
    "4" "Call Variants (Diploid)" \
    "5" "Filter Variants" \
    "6" "Generate Statistics"
```

### Workflow 2: Polyploid Sample (Integrated)

```
┌─────────────────────────────────────┐
│   Polyploid Workflow (3n, 4n, ...)  │
├─────────────────────────────────────┤
│ 1. Quality Control (01_qc.sh)      │
│ 2. Alignment (02_align.sh)         │
│ 3. BAM Processing (03_process_bam) │
│ 4. Variant Calling (04a_triploid)  │
│    ├─ Includes filtering            │
│    └─ Includes FreeBayes            │
│ 5. Statistics (06_stats)           │
└─────────────────────────────────────┘
```

**Note:** Step 5 (filtering) is SKIPPED - it's integrated into step 4!

---

## Complete YAD Example Script

Here's a complete example for Step 4a (Polyploid):

```bash
#!/bin/bash

# Step 4a: Polyploid Variant Calling with Clair3 + FreeBayes

# Get parameters from YAD dialog
RESULT=$(yad --title="Step 4a: Polyploid Variant Calling" \
    --form --width=700 --height=500 \
    --text="Integrated Clair3 + FreeBayes pipeline for polyploid genomes\n\nThis includes variant discovery, filtering, and re-genotyping." \
    --field="Input BAM (sorted):FL" "" \
    --field="Reference FASTA:FL" "" \
    --field="Output Directory:DIR" "" \
    --field="Clair3 Model:CB" "r1041_e82_400bps_sup_v500!r1041_e82_400bps_sup_v410!r941_prom_sup_g5014!r941_prom_hac_g360+g422" \
    --field="Threads:NUM" "16!1..64!1" \
    --field="Sample Ploidy:NUM" "3!2..8!1" \
    --field="Filter Preset:CB" "Standard (20/5)!Relaxed (15/3)!Strict (30/10)!Publication (40/15)!Custom" \
    --field="  Custom Min QUAL:NUM" "20!10..50!1" \
    --field="  Custom Min Depth:NUM" "5!3..20!1" \
    --button="Run:0" \
    --button="Cancel:1")

# Check if user cancelled
if [ $? -ne 0 ]; then
    exit 1
fi

# Parse results
IFS='|' read -r BAM REFERENCE OUTDIR MODEL THREADS PLOIDY PRESET CUSTOM_QUAL CUSTOM_DEPTH <<< "$RESULT"

# Validate inputs
if [ -z "$BAM" ] || [ -z "$REFERENCE" ] || [ -z "$OUTDIR" ]; then
    yad --error --text="Please fill in all required fields!"
    exit 1
fi

# Determine filter values based on preset
case "$PRESET" in
    "Standard (20/5)")
        MIN_QUAL=20
        MIN_DEPTH=5
        ;;
    "Relaxed (15/3)")
        MIN_QUAL=15
        MIN_DEPTH=3
        ;;
    "Strict (30/10)")
        MIN_QUAL=30
        MIN_DEPTH=10
        ;;
    "Publication (40/15)")
        MIN_QUAL=40
        MIN_DEPTH=15
        ;;
    "Custom")
        MIN_QUAL=$CUSTOM_QUAL
        MIN_DEPTH=$CUSTOM_DEPTH
        ;;
esac

# Show confirmation dialog
yad --question --width=500 \
    --text="Ready to run variant calling with these settings:\n\n\
BAM: $BAM\n\
Reference: $REFERENCE\n\
Output: $OUTDIR\n\
Model: $MODEL\n\
Ploidy: ${PLOIDY}n\n\
Filters: QUAL≥${MIN_QUAL}, DP≥${MIN_DEPTH}\n\
Threads: $THREADS\n\n\
This may take several hours. Continue?" \
    --button="Yes:0" --button="No:1"

if [ $? -ne 0 ]; then
    exit 1
fi

# Run the pipeline with progress bar
(
    echo "10" ; echo "# Running Clair3 variant discovery..."
    ./04a_call_variants_clair3_triploid_v2.sh \
        "$BAM" "$REFERENCE" "$OUTDIR" "$MODEL" "$THREADS" "$PLOIDY" "$MIN_QUAL" "$MIN_DEPTH"
    
    if [ $? -eq 0 ]; then
        echo "100" ; echo "# Complete!"
    else
        echo "100" ; echo "# ERROR: Pipeline failed"
        exit 1
    fi
) | yad --progress --title="Running Variant Calling" \
    --width=500 --height=100 \
    --auto-close --auto-kill \
    --no-cancel

# Check result
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    # Success - show results
    FINAL_VCF="${OUTDIR}/03_final_genotypes/final_ploidy${PLOIDY}.vcf.gz"
    SUMMARY="${OUTDIR}/pipeline_summary.txt"
    
    yad --text-info --title="Pipeline Complete" \
        --width=700 --height=500 \
        --filename="$SUMMARY" \
        --button="Open Output Folder:bash -c 'xdg-open ${OUTDIR}'" \
        --button="View VCF:bash -c 'bcftools view ${FINAL_VCF} | less'" \
        --button="Close:0"
else
    # Failed - show error
    LOG_FILE="${OUTDIR}/triploid_pipeline.log"
    yad --text-info --title="Pipeline Failed - Error Log" \
        --width=700 --height=500 \
        --filename="$LOG_FILE" \
        --button="Close:1"
fi
```

---

## GUI Design Tips

### 1. **Step Validation**
Before running each script, verify:
- BAM is sorted and indexed (for variant calling steps)
- Reference FASTA exists and is indexed
- Output directories are writable
- Sufficient disk space available

```bash
# Example validation
if [ ! -f "${BAM}.bai" ]; then
    yad --error --text="BAM file must be indexed!\n\nRun Step 3 (Process BAM) first."
    exit 1
fi
```

### 2. **Real-Time Progress**
Show script output in real-time:

```bash
./script.sh "$ARG1" "$ARG2" 2>&1 | \
    yad --text-info --title="Running..." \
    --width=800 --height=600 \
    --tail --button="Stop:1" --button="Background:0"
```

### 3. **File Browser Defaults**
Set sensible defaults for file browsers:

```bash
--field="Input BAM:FL" "$LAST_BAM_DIR" \
--field="Reference:FL" "$REFERENCE_GENOME"
```

Save user selections to a config file for next time.

### 4. **Tooltips**
Add helpful tooltips:

```bash
--field="Ploidy:NUM" "3!2..8!1" \
--field="":LBL "" \
--text="<b>Ploidy Guide:</b>\n\
2 = Diploid (most animals)\n\
3 = Triploid (banana, seedless watermelon)\n\
4 = Tetraploid (potato, cotton)\n\
6 = Hexaploid (wheat)"
```

### 5. **Pipeline State Tracking**
Track which steps have been completed:

```bash
# Save state
echo "BAM_SORTED=1" >> pipeline_state.txt
echo "VARIANTS_CALLED=1" >> pipeline_state.txt

# Load state and enable/disable buttons
source pipeline_state.txt
if [ "$BAM_SORTED" = "1" ]; then
    # Enable variant calling button
fi
```

### 6. **Multi-Sample Mode**
Add batch processing option:

```bash
# Sample list dialog
yad --list --title="Select Samples" \
    --column="Select:CHK" --column="Sample" --column="BAM" \
    TRUE "Sample1" "/path/sample1.bam" \
    TRUE "Sample2" "/path/sample2.bam" \
    FALSE "Sample3" "/path/sample3.bam"
```

---

## Error Handling

### Common Errors and User-Friendly Messages

```bash
# Example error handler
handle_error() {
    local exit_code=$1
    local log_file=$2
    
    case $exit_code in
        1)
            yad --error --text="Input file not found.\n\nCheck that all files exist and paths are correct."
            ;;
        2)
            yad --error --text="BAM file not indexed.\n\nPlease run Step 3 (Process BAM) first."
            ;;
        137)
            yad --error --text="Pipeline ran out of memory.\n\nTry reducing thread count or use a machine with more RAM."
            ;;
        *)
            yad --error --text="Pipeline failed with error code $exit_code.\n\nCheck log file:\n$log_file"
            ;;
    esac
}

# Usage
./script.sh "$ARGS" || handle_error $? "$LOG_FILE"
```

---

## Sample YAD Main Menu

```bash
#!/bin/bash

# ONT Variant Calling Pipeline - Main Menu

while true; do
    CHOICE=$(yad --list --title="ONT Variant Calling Pipeline" \
        --width=600 --height=500 \
        --text="<b>Oxford Nanopore Variant Calling Pipeline v2.0</b>\n\nSelect a step to run:" \
        --column="Step" --column="Name" --column="Status" \
        "1" "Quality Control" "✓" \
        "2" "Read Alignment" "✓" \
        "3" "Process BAM" "✓" \
        "4a" "Variant Calling (Diploid)" "○" \
        "4b" "Variant Calling (Polyploid)" "○" \
        "4c" "Variant Calling (Longshot)" "○" \
        "5" "Filter Variants" "○" \
        "6" "Variant Statistics" "○" \
        "---" "---" "---" \
        "CMP2" "Compare 2 VCFs" "○" \
        "CMP3" "Compare 3 VCFs" "○" \
        --button="Run Selected:0" \
        --button="Settings:2" \
        --button="Exit:1")
    
    ret=$?
    
    case $ret in
        1)
            # Exit
            break
            ;;
        2)
            # Settings
            show_settings_dialog
            ;;
        0)
            # Run selected step
            STEP=$(echo "$CHOICE" | cut -d'|' -f1)
            run_step "$STEP"
            ;;
    esac
done
```

---

## Configuration Management

Save and load user preferences:

```bash
CONFIG_FILE="$HOME/.config/ont_pipeline/config"

# Save settings
save_config() {
    mkdir -p "$(dirname $CONFIG_FILE)"
    cat > "$CONFIG_FILE" <<EOF
REFERENCE_GENOME="$REFERENCE"
DEFAULT_THREADS=$THREADS
DEFAULT_MODEL=$MODEL
LAST_OUTPUT_DIR=$OUTDIR
EOF
}

# Load settings
load_config() {
    if [ -f "$CONFIG_FILE" ]; then
        source "$CONFIG_FILE"
    fi
}
```

---

## Performance Tips

1. **Thread Count:** Default to physical cores, max at physical cores × 2
2. **Memory Warning:** Alert if <256GB RAM for large genomes
3. **Disk Space:** Check before running (need ~500GB per sample)
4. **Background Jobs:** Allow running multiple steps in background

---

## Final Notes

- All scripts are standalone - no dependencies between them except data files
- Scripts validate inputs before running
- All scripts write detailed logs
- Tab-separated output files are importable into Excel, R, Python
- Scripts handle both compressed (.gz) and uncompressed files

For questions or issues, check script log files first - they contain detailed error messages.

---

**Last Updated:** January 2, 2025
